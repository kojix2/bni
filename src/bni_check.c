#include "bni_internal.h"

#include <errno.h>
#include <getopt.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/tbx.h> /* hts_get_bgzfp() is declared here in htslib */

enum {
  OPT_QUICK = 1000,
  OPT_FULL = 1001,
};

static void usage_check(FILE *fp) {
  (void)fprintf(
      fp, "Usage:\n"
          "  bni check [options] <in.name.bam>\n\n"
          "Options:\n"
          "  -i, --index FILE       index file [default: <in.bam>.bni]\n"
          "      --quick            check metadata only [default]\n"
          "      --full             scan indexed BGZF-block ranges and verify boundaries/counts\n"
          "  -@, --threads INT      decompression threads\n"
          "  -h, --help             show this help\n");
}

static uint64_t header_hash64(sam_hdr_t *hdr) {
  const char *text = sam_hdr_str(hdr);
  size_t len = sam_hdr_length(hdr);
  if (text == NULL || len == SIZE_MAX) {
    return bni_fnv1a64("", 0);
  }
  return bni_fnv1a64(text, len);
}

static int check_metadata(const bni_index_t *idx, const char *bam_path, sam_hdr_t *hdr) {
  uint64_t bam_size = 0;
  int64_t bam_mtime = 0;
  int ok = 0;
  if (bni_file_metadata(bam_path, &bam_size, &bam_mtime) != 0) {
    bni_print_error("could not stat %s: %s", bam_path, strerror(errno));
    return -1;
  }
  if (idx->header.bam_size != bam_size) {
    bni_print_error("BAM size mismatch: index=%" PRIu64 " actual=%" PRIu64, idx->header.bam_size,
                    bam_size);
    ok = -1;
  }
  if (idx->header.bam_mtime != bam_mtime) {
    bni_print_error("BAM mtime mismatch: index=%" PRId64 " actual=%" PRId64, idx->header.bam_mtime,
                    bam_mtime);
    ok = -1;
  }
  uint64_t hh = header_hash64(hdr);
  if (idx->header.header_hash != hh) {
    bni_print_error("BAM header hash mismatch: index=%016" PRIx64 " actual=%016" PRIx64,
                    idx->header.header_hash, hh);
    ok = -1;
  }
  return ok;
}

typedef struct {
  samFile *in;
  sam_hdr_t *hdr;
  BGZF *bgzf_fp;
  bam1_t *record;
} full_check_ctx_t;

static int resolve_entry_names(const bni_index_t *idx, const bni_entry_t *entry,
                               uint64_t entry_index, const char **first_out,
                               const char **last_out) {
  const char *first = bni_entry_first_name(idx, entry);
  const char *last = bni_entry_last_name(idx, entry);
  if (first == NULL || last == NULL) {
    bni_print_error("entry %" PRIu64 " has invalid QNAME offsets", entry_index);
    return -1;
  }
  *first_out = first;
  *last_out = last;
  return 0;
}

static int check_record_position(const full_check_ctx_t *ctx, const bni_entry_t *entry,
                                 uint64_t entry_index, const char *first, const char *last) {
  int64_t pos = bgzf_tell(ctx->bgzf_fp);
  if (pos < 0) {
    bni_print_error("bgzf_tell failed for entry %" PRIu64 " (%s..%s)", entry_index, first, last);
    return -1;
  }
  if ((uint64_t)pos >= entry->end_voff) {
    bni_print_error("entry %" PRIu64 " ended before %u records", entry_index, entry->n_records);
    return -1;
  }
  return 0;
}

static int check_entry_start(const full_check_ctx_t *ctx, const bni_entry_t *entry,
                             uint64_t entry_index, const char *first, const char *last) {
  int64_t pos = bgzf_tell(ctx->bgzf_fp);
  if (pos < 0) {
    bni_print_error("bgzf_tell failed before entry %" PRIu64 " (%s..%s)", entry_index, first, last);
    return -1;
  }
  if ((uint64_t)pos != entry->beg_voff) {
    bni_print_error("entry %" PRIu64 " beg_voff mismatch: index=%" PRIu64 " BAM=%" PRIu64,
                    entry_index, entry->beg_voff, (uint64_t)pos);
    return -1;
  }
  return 0;
}

static int validate_record_name(const bni_entry_t *entry, uint64_t entry_index,
                                uint32_t record_index, const char *first, const char *last,
                                const char *previous, const char *qname) {
  if (record_index == 0 && strcmp(qname, first) != 0) {
    bni_print_error("entry %" PRIu64 " first QNAME mismatch: index='%s' BAM='%s'", entry_index,
                    first, qname);
    return -1;
  }
  if (record_index + 1 == entry->n_records && strcmp(qname, last) != 0) {
    bni_print_error("entry %" PRIu64 " last QNAME mismatch: index='%s' BAM='%s'", entry_index, last,
                    qname);
    return -1;
  }
  if (previous != NULL && strcmp(previous, qname) > 0) {
    bni_print_error("BAM QNAME order decreases inside entry %" PRIu64 " near '%s' -> '%s'",
                    entry_index, previous, qname);
    return -1;
  }
  return 0;
}

static int remember_qname(char **previous, const char *qname, uint64_t entry_index) {
  char *copy = bni_strdup(qname);
  if (copy == NULL) {
    bni_print_error("out of memory while checking entry %" PRIu64, entry_index);
    return -1;
  }
  free(*previous);
  *previous = copy;
  return 0;
}

static int check_entry_record(const full_check_ctx_t *ctx, const bni_entry_t *entry,
                              uint64_t entry_index, uint32_t record_index, const char *first,
                              const char *last, char **previous) {
  if (check_record_position(ctx, entry, entry_index, first, last) != 0) {
    return -1;
  }
  int ret = sam_read1(ctx->in, ctx->hdr, ctx->record);
  if (ret < 0) {
    bni_print_error("unexpected EOF while checking entry %" PRIu64 " (%s..%s)", entry_index, first,
                    last);
    return -1;
  }
  const char *qname = bam_get_qname(ctx->record);
  if (qname == NULL) {
    bni_print_error("entry %" PRIu64 " contains record with NULL QNAME", entry_index);
    return -1;
  }
  if (validate_record_name(entry, entry_index, record_index, first, last, *previous, qname) != 0) {
    return -1;
  }
  return remember_qname(previous, qname, entry_index);
}

static int check_entry_end(const full_check_ctx_t *ctx, const bni_entry_t *entry,
                           uint64_t entry_index) {
  int64_t pos_after = bgzf_tell(ctx->bgzf_fp);
  if (pos_after < 0) {
    bni_print_error("bgzf_tell failed after entry %" PRIu64, entry_index);
    return -1;
  }
  if ((uint64_t)pos_after != entry->end_voff) {
    bni_print_error("entry %" PRIu64 " end_voff mismatch: index=%" PRIu64 " BAM=%" PRIu64,
                    entry_index, entry->end_voff, (uint64_t)pos_after);
    return -1;
  }
  return 0;
}

static int check_index_entry(const full_check_ctx_t *ctx, const bni_index_t *idx,
                             uint64_t entry_index) {
  const bni_entry_t *entry = &idx->entries[entry_index];
  const char *first = NULL;
  const char *last = NULL;
  if (resolve_entry_names(idx, entry, entry_index, &first, &last) != 0) {
    return -1;
  }
  if (check_entry_start(ctx, entry, entry_index, first, last) != 0) {
    return -1;
  }

  char *previous = NULL;
  for (uint32_t record_index = 0; record_index < entry->n_records; ++record_index) {
    if (check_entry_record(ctx, entry, entry_index, record_index, first, last, &previous) != 0) {
      free(previous);
      return -1;
    }
  }
  free(previous);
  return check_entry_end(ctx, entry, entry_index);
}

static int seek_full_check_start(const full_check_ctx_t *ctx, const bni_index_t *idx) {
  if (idx->header.n_blocks == 0) {
    return 0;
  }
  const bni_entry_t *first_entry = &idx->entries[0];
  if (bgzf_seek(ctx->bgzf_fp, (int64_t)first_entry->beg_voff, SEEK_SET) < 0) {
    bni_print_error("bgzf_seek failed for first entry");
    return -1;
  }
  return 0;
}

static int full_check(const bni_index_t *idx, samFile *in, sam_hdr_t *hdr, BGZF *bgzf_fp) {
  full_check_ctx_t ctx = {.in = in, .hdr = hdr, .bgzf_fp = bgzf_fp, .record = bam_init1()};
  if (ctx.record == NULL) {
    bni_print_error("failed to allocate BAM record");
    return -1;
  }
  if (seek_full_check_start(&ctx, idx) != 0) {
    bam_destroy1(ctx.record);
    return -1;
  }
  for (uint64_t entry_index = 0; entry_index < idx->header.n_blocks; ++entry_index) {
    if (check_index_entry(&ctx, idx, entry_index) != 0) {
      bam_destroy1(ctx.record);
      return -1;
    }
  }
  bam_destroy1(ctx.record);
  return 0;
}

typedef struct {
  const char *index_path_arg;
  int do_full;
  int threads;
} check_options_t;

static int parse_check_threads_arg(int argc, char **argv, int *threads) {
  if (optarg == NULL) {
    if (optind < argc && bni_parse_threads(argv[optind], threads) == 0) {
      optind++;
      return 0;
    }
    bni_print_error("missing or invalid argument for -@/--threads");
    return -1;
  }
  if (bni_parse_threads(optarg, threads) != 0) {
    bni_print_error("invalid thread count '%s'", optarg);
    return -1;
  }
  return 0;
}

static int parse_check_options(int argc, char **argv, check_options_t *options) {
  static const struct option long_opts[] = {
      {"index", required_argument, NULL, 'i'}, {"quick", no_argument, NULL, OPT_QUICK},
      {"full", no_argument, NULL, OPT_FULL},   {"threads", required_argument, NULL, '@'},
      {"help", no_argument, NULL, 'h'},        {0, 0, 0, 0}};
  optind = 1;
  int c;
  while ((c = getopt_long(argc, argv, "i:@::h", long_opts, NULL)) != -1) {
    switch (c) {
    case 'i':
      options->index_path_arg = optarg;
      break;
    case '@':
      if (parse_check_threads_arg(argc, argv, &options->threads) != 0) {
        return -1;
      }
      break;
    case OPT_QUICK:
      options->do_full = 0;
      break;
    case OPT_FULL:
      options->do_full = 1;
      break;
    case 'h':
      usage_check(stdout);
      return 1;
    default:
      usage_check(stderr);
      return -1;
    }
  }
  return 0;
}

static int run_requested_check(const bni_index_t *idx, const char *bam_path, samFile *in,
                               sam_hdr_t *hdr, int do_full) {
  int status = 0;
  if (check_metadata(idx, bam_path, hdr) != 0) {
    status = 1;
  }
  if (status == 0 && do_full) {
    BGZF *bgzf_fp = hts_get_bgzfp(in);
    if (bgzf_fp == NULL) {
      bni_print_error("failed to access BGZF handle");
      status = 1;
    } else if (full_check(idx, in, hdr, bgzf_fp) != 0) {
      status = 1;
    }
  }
  return status;
}

static int run_check_command(const char *bam_path, const check_options_t *options) {
  char *default_index = NULL;
  const char *index_path = options->index_path_arg;
  if (index_path == NULL) {
    default_index = bni_default_index_path(bam_path);
    if (default_index == NULL) {
      bni_print_error("out of memory while building default index path");
      return 1;
    }
    index_path = default_index;
  }
  bni_index_t idx;
  if (bni_load_index_file(index_path, &idx) != 0) {
    free(default_index);
    return 1;
  }
  samFile *alignment = sam_open(bam_path, "r");
  if (alignment == NULL) {
    bni_print_error("could not open %s", bam_path);
    bni_index_destroy(&idx);
    free(default_index);
    return 1;
  }
  if (options->threads > 0 && hts_set_threads(alignment, options->threads) != 0) {
    bni_print_warning("failed to enable input threads");
  }
  sam_hdr_t *hdr = sam_hdr_read(alignment);
  if (hdr == NULL) {
    bni_print_error("failed to read BAM header from %s", bam_path);
    sam_close(alignment);
    bni_index_destroy(&idx);
    free(default_index);
    return 1;
  }
  int status = run_requested_check(&idx, bam_path, alignment, hdr, options->do_full);
  if (status == 0) {
    printf("OK\n");
  }
  sam_hdr_destroy(hdr);
  sam_close(alignment);
  bni_index_destroy(&idx);
  free(default_index);
  return status;
}

int bni_cmd_check(int argc, char **argv) {
  check_options_t options = {0};
  int parse_status = parse_check_options(argc, argv, &options);
  if (parse_status > 0) {
    return 0;
  }
  if (parse_status < 0) {
    return 1;
  }
  if (argc - optind != 1) {
    usage_check(stderr);
    return 1;
  }
  const char *bam_path = argv[optind];
  if (strcmp(bam_path, "-") == 0) {
    bni_print_error("checking stdin is not supported");
    return 1;
  }
  return run_check_command(bam_path, &options);
}
