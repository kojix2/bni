#define _POSIX_C_SOURCE 200809L

#include "bni.h"

#include <errno.h>
#include <getopt.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>  /* hts_get_bgzfp() is declared here in htslib */

typedef enum { BNI_OUT_BAM, BNI_OUT_SAM, BNI_OUT_CRAM } out_format_t;

static void usage_get(FILE *fp) {
    fprintf(fp,
        "Usage:\n"
        "  bni get [options] <in.name.bam> <read-name>\n"
        "  bni get [options] -f names.txt <in.name.bam>\n\n"
        "Extract all BAM records for read names using <in.name.bam>.bni.\n\n"
        "Options:\n"
        "  -i, --index FILE           index file [default: <in.bam>.bni]\n"
        "  -o, --output FILE          output file [default: stdout]\n"
        "  -O, --output-format FMT    bam|sam|cram [default: bam, or inferred from -o]\n"
        "  -f, --name-file FILE       read names, one per line\n"
        "  -@, --threads INT          input/output threads\n"
        "      --no-header            do not write SAM header; valid only with -O sam\n"
        "      --with-header          write header [default]\n"
        "      --missing-ok           exit 0 even if some names are absent\n"
        "      --list-missing         print missing names to stderr\n"
        "      --ignore-metadata      do not require BAM size/mtime/header hash to match index\n"
        "  -h, --help                show this help\n");
}

static int parse_output_format(const char *s, out_format_t *fmt) {
    if (strcmp(s, "bam") == 0) { *fmt = BNI_OUT_BAM; return 0; }
    if (strcmp(s, "sam") == 0) { *fmt = BNI_OUT_SAM; return 0; }
    if (strcmp(s, "cram") == 0) { *fmt = BNI_OUT_CRAM; return 0; }
    return -1;
}

static out_format_t infer_output_format(const char *path) {
    if (path != NULL) {
        if (bni_has_suffix(path, ".sam")) return BNI_OUT_SAM;
        if (bni_has_suffix(path, ".cram")) return BNI_OUT_CRAM;
    }
    return BNI_OUT_BAM;
}

static const char *mode_for_format(out_format_t fmt) {
    switch (fmt) { case BNI_OUT_BAM: return "wb"; case BNI_OUT_SAM: return "w"; case BNI_OUT_CRAM: return "wc"; }
    return "wb";
}

static uint64_t header_hash64(sam_hdr_t *hdr) {
    const char *text = sam_hdr_str(hdr);
    size_t len = sam_hdr_length(hdr);
    if (text == NULL || len == SIZE_MAX) return bni_fnv1a64("", 0);
    return bni_fnv1a64(text, len);
}

static int check_metadata(const bni_index_t *idx, const char *bam_path, sam_hdr_t *hdr) {
    uint64_t bam_size = 0; int64_t bam_mtime = 0;
    if (bni_file_metadata(bam_path, &bam_size, &bam_mtime) != 0) { bni_print_error("could not stat %s: %s", bam_path, strerror(errno)); return -1; }
    if (idx->header.bam_size != bam_size) { bni_print_error("BAM size differs from index metadata; rebuild the BNI index or use --ignore-metadata"); return -1; }
    if (idx->header.bam_mtime != bam_mtime) { bni_print_error("BAM mtime differs from index metadata; rebuild the BNI index or use --ignore-metadata"); return -1; }
    uint64_t hh = header_hash64(hdr);
    if (idx->header.header_hash != hh) { bni_print_error("BAM header hash differs from index metadata; rebuild the BNI index or use --ignore-metadata"); return -1; }
    return 0;
}

static char *trim_line(char *line) {
    if (line == NULL) return NULL;
    size_t n = strlen(line);
    while (n > 0 && (line[n - 1] == '\n' || line[n - 1] == '\r')) line[--n] = '\0';
    return line;
}

static int write_one_name(const bni_index_t *idx, samFile *in, sam_hdr_t *hdr, BGZF *bgzf_fp,
                          samFile *out, bam1_t *b, const char *name,
                          int list_missing, uint64_t *missing_count) {
    const bni_entry_t *entry = bni_find_entry(idx, name);
    if (entry == NULL) { if (list_missing) fprintf(stderr, "%s\n", name); (*missing_count)++; return 0; }
    if (bgzf_seek(bgzf_fp, (int64_t)entry->beg_voff, SEEK_SET) < 0) { bni_print_error("bgzf_seek failed for QNAME '%s'", name); return -1; }
    uint32_t seen = 0;
    while (1) {
        int64_t pos = bgzf_tell(bgzf_fp);
        if (pos < 0) { bni_print_error("bgzf_tell failed while reading QNAME '%s'", name); return -1; }
        if ((uint64_t)pos >= entry->end_voff) break;
        int ret = sam_read1(in, hdr, b);
        if (ret < 0) { bni_print_error("unexpected end of BAM while reading QNAME '%s'", name); return -1; }
        const char *qname = bam_get_qname(b);
        if (qname == NULL || strcmp(qname, name) != 0) { bni_print_error("index range for '%s' contains record with QNAME '%s'; rebuild the index", name, qname ? qname : "(null)"); return -1; }
        if (sam_write1(out, hdr, b) < 0) { bni_print_error("failed writing output record for QNAME '%s'", name); return -1; }
        seen++;
    }
    if (seen != entry->n_records) { bni_print_error("index range for '%s' contained %u records; expected %u", name, seen, entry->n_records); return -1; }
    return 0;
}

int bni_cmd_get(int argc, char **argv) {
    const char *index_path_arg = NULL, *out_path = NULL, *names_path = NULL, *fmt_arg = NULL;
    int write_header = 1, missing_ok = 0, list_missing = 0, ignore_metadata = 0, threads = 0;
    static const struct option long_opts[] = {
        {"index", required_argument, NULL, 'i'}, {"output", required_argument, NULL, 'o'},
        {"output-format", required_argument, NULL, 'O'}, {"name-file", required_argument, NULL, 'f'},
        {"threads", required_argument, NULL, '@'}, {"no-header", no_argument, NULL, 1000},
        {"with-header", no_argument, NULL, 1001}, {"missing-ok", no_argument, NULL, 1002},
        {"list-missing", no_argument, NULL, 1003}, {"ignore-metadata", no_argument, NULL, 1004},
        {"help", no_argument, NULL, 'h'}, {0,0,0,0}
    };
    optind = 1;
    int c;
    while ((c = getopt_long(argc, argv, "i:o:O:f:@::h", long_opts, NULL)) != -1) {
        switch (c) {
        case 'i': index_path_arg = optarg; break;
        case 'o': out_path = optarg; break;
        case 'O': fmt_arg = optarg; break;
        case 'f': names_path = optarg; break;
        case '@':
            if (optarg == NULL) {
                if (optind < argc && bni_parse_threads(argv[optind], &threads) == 0) { optind++; break; }
                bni_print_error("missing or invalid argument for -@/--threads"); return 1;
            }
            if (bni_parse_threads(optarg, &threads) != 0) { bni_print_error("invalid thread count '%s'", optarg); return 1; }
            break;
        case 1000: write_header = 0; break;
        case 1001: write_header = 1; break;
        case 1002: missing_ok = 1; break;
        case 1003: list_missing = 1; break;
        case 1004: ignore_metadata = 1; break;
        case 'h': usage_get(stdout); return 0;
        default: usage_get(stderr); return 1;
        }
    }
    const char *bam_path = NULL, *single_name = NULL;
    if (names_path != NULL) { if (argc - optind != 1) { usage_get(stderr); return 1; } bam_path = argv[optind]; }
    else { if (argc - optind != 2) { usage_get(stderr); return 1; } bam_path = argv[optind]; single_name = argv[optind + 1]; }
    if (strcmp(bam_path, "-") == 0) { bni_print_error("random access from stdin is not supported"); return 1; }
    out_format_t out_fmt = infer_output_format(out_path);
    if (fmt_arg != NULL && parse_output_format(fmt_arg, &out_fmt) != 0) { bni_print_error("unknown output format '%s'; expected bam, sam, or cram", fmt_arg); return 1; }
    if (!write_header && out_fmt != BNI_OUT_SAM) { bni_print_error("--no-header is only supported for SAM output"); return 1; }
    char *default_index = NULL;
    const char *index_path = index_path_arg;
    if (index_path == NULL) { default_index = bni_default_index_path(bam_path); if (default_index == NULL) { bni_print_error("out of memory while building default index path"); return 1; } index_path = default_index; }
    bni_index_t idx;
    if (bni_load_index_file(index_path, &idx) != 0) { free(default_index); return 1; }
    samFile *in = sam_open(bam_path, "r");
    if (in == NULL) { bni_print_error("could not open %s", bam_path); bni_index_destroy(&idx); free(default_index); return 1; }
    if (threads > 0 && hts_set_threads(in, threads) != 0) bni_print_warning("failed to enable input threads");
    sam_hdr_t *hdr = sam_hdr_read(in);
    if (hdr == NULL) { bni_print_error("failed to read BAM header from %s", bam_path); sam_close(in); bni_index_destroy(&idx); free(default_index); return 1; }
    const htsFormat *fmt = hts_get_format(in);
    if (fmt != NULL && fmt->format != bam) { bni_print_error("input is not BAM; BNI v1 supports BGZF-compressed BAM only"); sam_hdr_destroy(hdr); sam_close(in); bni_index_destroy(&idx); free(default_index); return 1; }
    if (!ignore_metadata && check_metadata(&idx, bam_path, hdr) != 0) { sam_hdr_destroy(hdr); sam_close(in); bni_index_destroy(&idx); free(default_index); return 1; }
    BGZF *bgzf_fp = hts_get_bgzfp(in);
    if (bgzf_fp == NULL) { bni_print_error("failed to access BGZF handle"); sam_hdr_destroy(hdr); sam_close(in); bni_index_destroy(&idx); free(default_index); return 1; }
    samFile *out = sam_open(out_path ? out_path : "-", mode_for_format(out_fmt));
    if (out == NULL) { bni_print_error("could not open output"); sam_hdr_destroy(hdr); sam_close(in); bni_index_destroy(&idx); free(default_index); return 1; }
    if (threads > 0 && hts_set_threads(out, threads) != 0) bni_print_warning("failed to enable output threads");
    if (write_header && sam_hdr_write(out, hdr) != 0) { bni_print_error("failed to write output header"); sam_close(out); sam_hdr_destroy(hdr); sam_close(in); bni_index_destroy(&idx); free(default_index); return 1; }
    bam1_t *b = bam_init1();
    if (b == NULL) { bni_print_error("failed to allocate BAM record"); sam_close(out); sam_hdr_destroy(hdr); sam_close(in); bni_index_destroy(&idx); free(default_index); return 1; }
    int status = 0; uint64_t missing = 0;
    if (names_path != NULL) {
        FILE *nf = strcmp(names_path, "-") == 0 ? stdin : fopen(names_path, "r");
        if (nf == NULL) { bni_print_error("could not open name file %s: %s", names_path, strerror(errno)); status = 1; }
        else {
            char *line = NULL; size_t cap = 0; ssize_t nread;
            while ((nread = getline(&line, &cap, nf)) >= 0) { (void)nread; char *name = trim_line(line); if (name[0] == '\0') continue; if (write_one_name(&idx, in, hdr, bgzf_fp, out, b, name, list_missing, &missing) != 0) { status = 1; break; } }
            free(line); if (nf != stdin) fclose(nf);
        }
    } else {
        if (write_one_name(&idx, in, hdr, bgzf_fp, out, b, single_name, list_missing, &missing) != 0) status = 1;
    }
    if (sam_close(out) != 0) { bni_print_error("failed closing output"); status = 1; }
    bam_destroy1(b); sam_hdr_destroy(hdr); sam_close(in); bni_index_destroy(&idx); free(default_index);
    if (status == 0 && missing > 0 && !missing_ok) { char mbuf[64]; bni_format_u64(mbuf, sizeof(mbuf), missing); bni_print_error("%s requested read name(s) were not found; use --missing-ok to ignore", mbuf); status = 1; }
    return status;
}
