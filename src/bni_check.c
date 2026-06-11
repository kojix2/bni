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
#include <htslib/sam.h>
#include <htslib/tbx.h>  /* hts_get_bgzfp() is declared here in htslib */

static void usage_check(FILE *fp) {
    fprintf(fp,
        "Usage:\n"
        "  bni check [options] <in.name.bam>\n\n"
        "Options:\n"
        "  -i, --index FILE       index file [default: <in.bam>.bni]\n"
        "      --quick            check metadata only [default]\n"
        "      --full             also seek every indexed range and verify QNAME/count\n"
        "  -@, --threads INT      decompression threads\n"
        "  -h, --help             show this help\n");
}

static uint64_t header_hash64(sam_hdr_t *hdr) {
    const char *text = sam_hdr_str(hdr);
    size_t len = sam_hdr_length(hdr);
    if (text == NULL || len == SIZE_MAX) return bni_fnv1a64("", 0);
    return bni_fnv1a64(text, len);
}

static int check_metadata(const bni_index_t *idx, const char *bam_path, sam_hdr_t *hdr) {
    uint64_t bam_size = 0; int64_t bam_mtime = 0; int ok = 0;
    if (bni_file_metadata(bam_path, &bam_size, &bam_mtime) != 0) { bni_print_error("could not stat %s: %s", bam_path, strerror(errno)); return -1; }
    if (idx->header.bam_size != bam_size) { bni_print_error("BAM size mismatch: index=%" PRIu64 " actual=%" PRIu64, idx->header.bam_size, bam_size); ok = -1; }
    if (idx->header.bam_mtime != bam_mtime) { bni_print_error("BAM mtime mismatch: index=%" PRId64 " actual=%" PRId64, idx->header.bam_mtime, bam_mtime); ok = -1; }
    uint64_t hh = header_hash64(hdr);
    if (idx->header.header_hash != hh) { bni_print_error("BAM header hash mismatch: index=%016" PRIx64 " actual=%016" PRIx64, idx->header.header_hash, hh); ok = -1; }
    return ok;
}

static int full_check(const bni_index_t *idx, samFile *in, sam_hdr_t *hdr, BGZF *bgzf_fp) {
    bam1_t *b = bam_init1();
    if (b == NULL) { bni_print_error("failed to allocate BAM record"); return -1; }
    for (uint64_t i = 0; i < idx->header.n_names; ++i) {
        const bni_entry_t *entry = &idx->entries[i];
        const char *name = bni_entry_name(idx, entry);
        if (name == NULL) { bni_print_error("entry %" PRIu64 " has invalid name", i); bam_destroy1(b); return -1; }
        if (bgzf_seek(bgzf_fp, (int64_t)entry->beg_voff, SEEK_SET) < 0) { bni_print_error("bgzf_seek failed for entry %" PRIu64 " (%s)", i, name); bam_destroy1(b); return -1; }
        uint32_t seen = 0;
        while (1) {
            int64_t pos = bgzf_tell(bgzf_fp);
            if (pos < 0) { bni_print_error("bgzf_tell failed for entry %" PRIu64 " (%s)", i, name); bam_destroy1(b); return -1; }
            if ((uint64_t)pos >= entry->end_voff) break;
            int ret = sam_read1(in, hdr, b);
            if (ret < 0) { bni_print_error("unexpected EOF while checking entry %" PRIu64 " (%s)", i, name); bam_destroy1(b); return -1; }
            const char *qname = bam_get_qname(b);
            if (qname == NULL || strcmp(qname, name) != 0) { bni_print_error("entry %" PRIu64 " expected QNAME '%s', saw '%s'", i, name, qname ? qname : "(null)"); bam_destroy1(b); return -1; }
            seen++;
        }
        if (seen != entry->n_records) { bni_print_error("entry %" PRIu64 " (%s) expected %u records, saw %u", i, name, entry->n_records, seen); bam_destroy1(b); return -1; }
    }
    bam_destroy1(b); return 0;
}

int bni_cmd_check(int argc, char **argv) {
    const char *index_path_arg = NULL; int do_full = 0, threads = 0;
    static const struct option long_opts[] = { {"index", required_argument, NULL, 'i'}, {"quick", no_argument, NULL, 1000}, {"full", no_argument, NULL, 1001}, {"threads", required_argument, NULL, '@'}, {"help", no_argument, NULL, 'h'}, {0,0,0,0} };
    optind = 1; int c;
    while ((c = getopt_long(argc, argv, "i:@::h", long_opts, NULL)) != -1) {
        switch (c) {
        case 'i': index_path_arg = optarg; break;
        case '@':
            if (optarg == NULL) {
                if (optind < argc && bni_parse_threads(argv[optind], &threads) == 0) { optind++; break; }
                bni_print_error("missing or invalid argument for -@/--threads"); return 1;
            }
            if (bni_parse_threads(optarg, &threads) != 0) { bni_print_error("invalid thread count '%s'", optarg); return 1; }
            break;
        case 1000: do_full = 0; break; case 1001: do_full = 1; break; case 'h': usage_check(stdout); return 0; default: usage_check(stderr); return 1;
        }
    }
    if (argc - optind != 1) { usage_check(stderr); return 1; }
    const char *bam_path = argv[optind];
    if (strcmp(bam_path, "-") == 0) { bni_print_error("checking stdin is not supported"); return 1; }
    char *default_index = NULL; const char *index_path = index_path_arg;
    if (index_path == NULL) { default_index = bni_default_index_path(bam_path); if (default_index == NULL) { bni_print_error("out of memory while building default index path"); return 1; } index_path = default_index; }
    bni_index_t idx;
    if (bni_load_index_file(index_path, &idx) != 0) { free(default_index); return 1; }
    samFile *in = sam_open(bam_path, "r");
    if (in == NULL) { bni_print_error("could not open %s", bam_path); bni_index_destroy(&idx); free(default_index); return 1; }
    if (threads > 0 && hts_set_threads(in, threads) != 0) bni_print_warning("failed to enable input threads");
    sam_hdr_t *hdr = sam_hdr_read(in);
    if (hdr == NULL) { bni_print_error("failed to read BAM header from %s", bam_path); sam_close(in); bni_index_destroy(&idx); free(default_index); return 1; }
    int status = 0;
    if (check_metadata(&idx, bam_path, hdr) != 0) status = 1;
    if (status == 0 && do_full) {
        BGZF *bgzf_fp = hts_get_bgzfp(in);
        if (bgzf_fp == NULL) { bni_print_error("failed to access BGZF handle"); status = 1; }
        else if (full_check(&idx, in, hdr, bgzf_fp) != 0) status = 1;
    }
    if (status == 0) printf("OK\n");
    sam_hdr_destroy(hdr); sam_close(in); bni_index_destroy(&idx); free(default_index); return status;
}
