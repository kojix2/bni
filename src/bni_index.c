#define _POSIX_C_SOURCE 200809L

#include "bni.h"

#include <errno.h>
#include <getopt.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>  /* hts_get_bgzfp() is declared here in htslib */

typedef struct { bni_entry_t *data; size_t len; size_t cap; } entry_vec_t;
typedef struct { char *data; size_t len; size_t cap; } strbuf_t;

static const char *g_sort_strings = NULL;

static void usage_index(FILE *fp) {
    fprintf(fp,
        "Usage:\n"
        "  bni index [options] <in.name.bam>\n\n"
        "Create <in.name.bam>.bni for a BAM sorted with samtools sort -N.\n\n"
        "Options:\n"
        "  -o, --output FILE          output index file [default: <in.bam>.bni]\n"
        "  -f, --force                overwrite an existing index\n"
        "  -@, --threads INT          decompression threads\n"
        "      --no-header-check      do not require @HD SO/SS tags\n"
        "  -v, --verbose              print progress summary\n"
        "  -h, --help                 show this help\n");
}

static void entry_vec_free(entry_vec_t *v) { free(v->data); v->data = NULL; v->len = v->cap = 0; }
static void strbuf_free(strbuf_t *s) { free(s->data); s->data = NULL; s->len = s->cap = 0; }

static int entry_vec_push(entry_vec_t *v, bni_entry_t e) {
    if (v->len == v->cap) {
        size_t new_cap = v->cap ? v->cap * 2 : 4096;
        if (new_cap < v->cap || new_cap > SIZE_MAX / sizeof(bni_entry_t)) { bni_print_error("too many read names for this platform"); return -1; }
        bni_entry_t *p = (bni_entry_t *)realloc(v->data, new_cap * sizeof(bni_entry_t));
        if (p == NULL) { bni_print_error("out of memory while growing entry table"); return -1; }
        v->data = p;
        v->cap = new_cap;
    }
    v->data[v->len++] = e;
    return 0;
}

static int strbuf_append_cstr(strbuf_t *s, const char *name, uint64_t *offset_out) {
    size_t n = strlen(name) + 1;
    if (s->len > UINT64_MAX) return -1;
    if (offset_out) *offset_out = (uint64_t)s->len;
    if (n > SIZE_MAX - s->len) { bni_print_error("string table is too large"); return -1; }
    size_t need = s->len + n;
    if (need > s->cap) {
        size_t new_cap = s->cap ? s->cap * 2 : 65536;
        while (new_cap < need) {
            if (new_cap > SIZE_MAX / 2) { new_cap = need; break; }
            new_cap *= 2;
        }
        char *p = (char *)realloc(s->data, new_cap);
        if (p == NULL) { bni_print_error("out of memory while growing string table"); return -1; }
        s->data = p;
        s->cap = new_cap;
    }
    memcpy(s->data + s->len, name, n);
    s->len += n;
    return 0;
}

static int add_group(entry_vec_t *entries, strbuf_t *strings, const char *name,
                     uint64_t beg_voff, uint64_t end_voff, uint32_t n_records) {
    uint64_t name_offset = 0;
    if (strbuf_append_cstr(strings, name, &name_offset) != 0) return -1;
    bni_entry_t e;
    memset(&e, 0, sizeof(e));
    e.name_offset = name_offset;
    e.beg_voff = beg_voff;
    e.end_voff = end_voff;
    e.n_records = n_records;
    return entry_vec_push(entries, e);
}

static int entry_cmp_by_name(const void *a, const void *b) {
    const bni_entry_t *ea = (const bni_entry_t *)a;
    const bni_entry_t *eb = (const bni_entry_t *)b;
    return strcmp(g_sort_strings + ea->name_offset, g_sort_strings + eb->name_offset);
}

static int header_is_queryname_lex(sam_hdr_t *hdr, int no_header_check) {
    if (no_header_check) return 0;
    kstring_t so = {0, 0, NULL};
    kstring_t ss = {0, 0, NULL};
    int ok = 0;
    if (sam_hdr_find_tag_hd(hdr, "SO", &so) != 0 || so.s == NULL || strcmp(so.s, "queryname") != 0) {
        bni_print_error("input header is not @HD SO:queryname; run: samtools sort -N -o out.name.bam in.bam");
        ok = -1; goto done;
    }
    if (sam_hdr_find_tag_hd(hdr, "SS", &ss) == 0 && ss.s != NULL) {
        if (strcmp(ss.s, "queryname:lexicographical") != 0) {
            bni_print_error("input header has @HD SS:%s, not SS:queryname:lexicographical", ss.s);
            ok = -1; goto done;
        }
    } else {
        bni_print_warning("@HD SS tag is absent; validating lexicographic QNAME order while indexing");
    }
done:
    free(so.s);
    free(ss.s);
    return ok;
}

static uint64_t header_hash64(sam_hdr_t *hdr) {
    const char *text = sam_hdr_str(hdr);
    size_t len = sam_hdr_length(hdr);
    if (text == NULL || len == SIZE_MAX) return bni_fnv1a64("", 0);
    return bni_fnv1a64(text, len);
}

int bni_cmd_index(int argc, char **argv) {
    const char *out_path_arg = NULL;
    int force = 0, threads = 0, no_header_check = 0, verbose = 0;
    static const struct option long_opts[] = {
        {"output", required_argument, NULL, 'o'}, {"force", no_argument, NULL, 'f'},
        {"threads", required_argument, NULL, '@'}, {"no-header-check", no_argument, NULL, 1000},
        {"verbose", no_argument, NULL, 'v'}, {"help", no_argument, NULL, 'h'}, {0,0,0,0}
    };
    optind = 1;
    int c;
    while ((c = getopt_long(argc, argv, "o:f@::vh", long_opts, NULL)) != -1) {
        switch (c) {
        case 'o': out_path_arg = optarg; break;
        case 'f': force = 1; break;
        case '@':
            if (optarg == NULL) {
                if (optind < argc && bni_parse_threads(argv[optind], &threads) == 0) { optind++; break; }
                bni_print_error("missing or invalid argument for -@/--threads"); return 1;
            }
            if (bni_parse_threads(optarg, &threads) != 0) { bni_print_error("invalid thread count '%s'", optarg); return 1; }
            break;
        case 'v': verbose = 1; break;
        case 1000: no_header_check = 1; break;
        case 'h': usage_index(stdout); return 0;
        default: usage_index(stderr); return 1;
        }
    }
    if (argc - optind != 1) { usage_index(stderr); return 1; }
    const char *bam_path = argv[optind];
    if (strcmp(bam_path, "-") == 0) { bni_print_error("indexing from stdin is not supported because BNI stores BGZF virtual offsets"); return 1; }
    char *default_out = NULL;
    const char *out_path = out_path_arg;
    if (out_path == NULL) {
        default_out = bni_default_index_path(bam_path);
        if (default_out == NULL) { bni_print_error("out of memory while building output path"); return 1; }
        out_path = default_out;
    }
    if (!force && bni_path_exists(out_path)) { bni_print_error("%s already exists; use --force to overwrite", out_path); free(default_out); return 1; }
    uint64_t bam_size = 0; int64_t bam_mtime = 0;
    if (bni_file_metadata(bam_path, &bam_size, &bam_mtime) != 0) { bni_print_error("could not stat %s: %s", bam_path, strerror(errno)); free(default_out); return 1; }

    samFile *fp = sam_open(bam_path, "r");
    if (fp == NULL) { bni_print_error("could not open %s", bam_path); free(default_out); return 1; }
    if (threads > 0 && hts_set_threads(fp, threads) != 0) bni_print_warning("failed to enable input threads");
    bam_hdr_t *hdr = sam_hdr_read(fp);
    if (hdr == NULL) { bni_print_error("failed to read BAM header from %s", bam_path); sam_close(fp); free(default_out); return 1; }
    const htsFormat *fmt = hts_get_format(fp);
    if (fmt != NULL && fmt->format != bam) { bni_print_error("input is not BAM; BNI v1 supports BGZF-compressed BAM only"); sam_hdr_destroy(hdr); sam_close(fp); free(default_out); return 1; }
    if (fmt != NULL && fmt->compression != bgzf) { bni_print_error("input is not BGZF-compressed BAM; BNI v1 cannot seek it safely"); sam_hdr_destroy(hdr); sam_close(fp); free(default_out); return 1; }
    if (header_is_queryname_lex(hdr, no_header_check) != 0) { sam_hdr_destroy(hdr); sam_close(fp); free(default_out); return 1; }
    BGZF *bgzf_fp = hts_get_bgzfp(fp);
    if (bgzf_fp == NULL) { bni_print_error("failed to access BGZF handle"); sam_hdr_destroy(hdr); sam_close(fp); free(default_out); return 1; }

    entry_vec_t entries = {0,0,0}; strbuf_t strings = {0,0,0};
    bam1_t *b = bam_init1();
    if (b == NULL) { bni_print_error("failed to allocate BAM record"); sam_hdr_destroy(hdr); sam_close(fp); free(default_out); return 1; }
    char *prev_name = NULL;
    uint64_t group_beg = 0, total_records = 0, last_rec_end = 0;
    uint32_t group_n = 0;
    int exit_status = 1;

    for (;;) {
        int64_t tell_before = bgzf_tell(bgzf_fp);
        if (tell_before < 0) { bni_print_error("bgzf_tell failed before reading record"); goto cleanup; }
        uint64_t rec_beg = (uint64_t)tell_before;
        int ret = sam_read1(fp, hdr, b);
        if (ret < 0) { if (ret < -1) { bni_print_error("error while reading BAM records"); goto cleanup; } break; }
        int64_t tell_after = bgzf_tell(bgzf_fp);
        if (tell_after < 0) { bni_print_error("bgzf_tell failed after reading record"); goto cleanup; }
        last_rec_end = (uint64_t)tell_after;
        const char *qname = bam_get_qname(b);
        if (qname == NULL || qname[0] == '\0') { bni_print_error("encountered record with empty QNAME"); goto cleanup; }
        if (prev_name == NULL) {
            prev_name = bni_strdup(qname);
            if (prev_name == NULL) { bni_print_error("out of memory while storing QNAME"); goto cleanup; }
            group_beg = rec_beg; group_n = 1;
        } else {
            int cmp = strcmp(prev_name, qname);
            if (cmp > 0) { bni_print_error("BAM is not lexicographically queryname-sorted near '%s' -> '%s'; use samtools sort -N", prev_name, qname); goto cleanup; }
            if (cmp == 0) {
                if (group_n == UINT32_MAX) { bni_print_error("too many records for QNAME '%s'", prev_name); goto cleanup; }
                group_n++;
            } else {
                if (add_group(&entries, &strings, prev_name, group_beg, rec_beg, group_n) != 0) goto cleanup;
                free(prev_name);
                prev_name = bni_strdup(qname);
                if (prev_name == NULL) { bni_print_error("out of memory while storing QNAME"); goto cleanup; }
                group_beg = rec_beg; group_n = 1;
            }
        }
        total_records++;
    }
    if (prev_name != NULL) {
        if (add_group(&entries, &strings, prev_name, group_beg, last_rec_end, group_n) != 0) goto cleanup;
    }
    if (entries.len > 1) {
        g_sort_strings = strings.data;
        qsort(entries.data, entries.len, sizeof(entries.data[0]), entry_cmp_by_name);
        g_sort_strings = NULL;
        for (size_t i = 1; i < entries.len; ++i) {
            const char *a = strings.data + entries.data[i - 1].name_offset;
            const char *bn = strings.data + entries.data[i].name_offset;
            if (strcmp(a, bn) == 0) { bni_print_error("QNAME '%s' appears in multiple non-contiguous groups", a); goto cleanup; }
        }
    }
    bni_file_header_t header; memset(&header, 0, sizeof(header));
    header.version = BNI_FORMAT_VERSION;
    header.header_size = BNI_HEADER_SIZE;
    header.flags = BNI_FLAG_QUERYNAME_GROUPED;
    header.n_names = (uint64_t)entries.len;
    header.n_records = total_records;
    header.entries_offset = BNI_HEADER_SIZE;
    header.strings_offset = BNI_HEADER_SIZE + (uint64_t)entries.len * (uint64_t)BNI_ENTRY_SIZE;
    header.strings_size = (uint64_t)strings.len;
    header.bam_size = bam_size;
    header.bam_mtime = bam_mtime;
    header.header_hash = header_hash64(hdr);
    header.sort_order = BNI_SORT_QUERYNAME_LEX;
    header.entry_size = BNI_ENTRY_SIZE;
    if (bni_write_index_file(out_path, &header, entries.data, strings.data ? strings.data : "") != 0) goto cleanup;
    if (verbose) {
        char nbuf[64], rbuf[64], sbuf[64];
        bni_format_u64(nbuf, sizeof(nbuf), header.n_names);
        bni_format_u64(rbuf, sizeof(rbuf), header.n_records);
        bni_format_u64(sbuf, sizeof(sbuf), header.strings_size);
        fprintf(stderr, "bni: wrote %s\n", out_path);
        fprintf(stderr, "bni: names=%s records=%s string_bytes=%s\n", nbuf, rbuf, sbuf);
    }
    exit_status = 0;
cleanup:
    free(prev_name); bam_destroy1(b); entry_vec_free(&entries); strbuf_free(&strings);
    sam_hdr_destroy(hdr); sam_close(fp); free(default_out);
    return exit_status;
}
