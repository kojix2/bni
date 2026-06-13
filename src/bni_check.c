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
      fp,
      "Usage:\n"
      "  bni check [options] <in.name.bam>\n\n"
      "Options:\n"
      "  -i, --index FILE       index file [default: <in.bam>.bni]\n"
      "      --quick            check metadata only [default]\n"
      "      --full             seek every indexed BGZF-block range and verify boundaries/counts\n"
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

static int full_check(const bni_index_t *idx, samFile *in, sam_hdr_t *hdr, BGZF *bgzf_fp) {
  bam1_t *b = bam_init1();
  if (b == NULL) {
    bni_print_error("failed to allocate BAM record");
    return -1;
  }
  for (uint64_t i = 0; i < idx->header.n_blocks; ++i) {
    const bni_entry_t *entry = &idx->entries[i];
    const char *first = bni_entry_first_name(idx, entry);
    const char *last = bni_entry_last_name(idx, entry);
    if (first == NULL || last == NULL) {
      bni_print_error("entry %" PRIu64 " has invalid QNAME offsets", i);
      bam_destroy1(b);
      return -1;
    }
    if (bgzf_seek(bgzf_fp, (int64_t)entry->beg_voff, SEEK_SET) < 0) {
      bni_print_error("bgzf_seek failed for entry %" PRIu64 " (%s..%s)", i, first, last);
      bam_destroy1(b);
      return -1;
    }
    char *prev = NULL;
    for (uint32_t j = 0; j < entry->n_records; ++j) {
      int64_t pos = bgzf_tell(bgzf_fp);
      if (pos < 0) {
        bni_print_error("bgzf_tell failed for entry %" PRIu64 " (%s..%s)", i, first, last);
        free(prev);
        bam_destroy1(b);
        return -1;
      }
      if ((uint64_t)pos >= entry->end_voff) {
        bni_print_error("entry %" PRIu64 " ended before %u records", i, entry->n_records);
        free(prev);
        bam_destroy1(b);
        return -1;
      }
      int ret = sam_read1(in, hdr, b);
      if (ret < 0) {
        bni_print_error("unexpected EOF while checking entry %" PRIu64 " (%s..%s)", i, first, last);
        free(prev);
        bam_destroy1(b);
        return -1;
      }
      const char *qname = bam_get_qname(b);
      if (qname == NULL) {
        bni_print_error("entry %" PRIu64 " contains record with NULL QNAME", i);
        free(prev);
        bam_destroy1(b);
        return -1;
      }
      if (j == 0 && strcmp(qname, first) != 0) {
        bni_print_error("entry %" PRIu64 " first QNAME mismatch: index='%s' BAM='%s'", i, first,
                        qname);
        free(prev);
        bam_destroy1(b);
        return -1;
      }
      if (j + 1 == entry->n_records && strcmp(qname, last) != 0) {
        bni_print_error("entry %" PRIu64 " last QNAME mismatch: index='%s' BAM='%s'", i, last,
                        qname);
        free(prev);
        bam_destroy1(b);
        return -1;
      }
      if (prev != NULL && strcmp(prev, qname) > 0) {
        bni_print_error("BAM QNAME order decreases inside entry %" PRIu64 " near '%s' -> '%s'", i,
                        prev, qname);
        free(prev);
        bam_destroy1(b);
        return -1;
      }
      free(prev);
      prev = bni_strdup(qname);
      if (prev == NULL) {
        bni_print_error("out of memory while checking entry %" PRIu64, i);
        bam_destroy1(b);
        return -1;
      }
    }
    free(prev);
    int64_t pos_after = bgzf_tell(bgzf_fp);
    if (pos_after < 0) {
      bni_print_error("bgzf_tell failed after entry %" PRIu64, i);
      bam_destroy1(b);
      return -1;
    }
    if ((uint64_t)pos_after != entry->end_voff) {
      bni_print_error("entry %" PRIu64 " end_voff mismatch: index=%" PRIu64 " BAM=%" PRIu64, i,
                      entry->end_voff, (uint64_t)pos_after);
      bam_destroy1(b);
      return -1;
    }
  }
  bam_destroy1(b);
  return 0;
}

int bni_cmd_check(int argc, char **argv) {
  const char *index_path_arg = NULL;
  int do_full = 0;
  int threads = 0;
  static const struct option long_opts[] = {
      {"index", required_argument, NULL, 'i'}, {"quick", no_argument, NULL, OPT_QUICK},
      {"full", no_argument, NULL, OPT_FULL},   {"threads", required_argument, NULL, '@'},
      {"help", no_argument, NULL, 'h'},        {0, 0, 0, 0}};
  optind = 1;
  int c;
  while ((c = getopt_long(argc, argv, "i:@::h", long_opts, NULL)) != -1) {
    switch (c) {
    case 'i':
      index_path_arg = optarg;
      break;
    case '@':
      if (optarg == NULL) {
        if (optind < argc && bni_parse_threads(argv[optind], &threads) == 0) {
          optind++;
          break;
        }
        bni_print_error("missing or invalid argument for -@/--threads");
        return 1;
      }
      if (bni_parse_threads(optarg, &threads) != 0) {
        bni_print_error("invalid thread count '%s'", optarg);
        return 1;
      }
      break;
    case OPT_QUICK:
      do_full = 0;
      break;
    case OPT_FULL:
      do_full = 1;
      break;
    case 'h':
      usage_check(stdout);
      return 0;
    default:
      usage_check(stderr);
      return 1;
    }
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
  char *default_index = NULL;
  const char *index_path = index_path_arg;
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
  samFile *in = sam_open(bam_path, "r");
  if (in == NULL) {
    bni_print_error("could not open %s", bam_path);
    bni_index_destroy(&idx);
    free(default_index);
    return 1;
  }
  if (threads > 0 && hts_set_threads(in, threads) != 0) {
    bni_print_warning("failed to enable input threads");
  }
  sam_hdr_t *hdr = sam_hdr_read(in);
  if (hdr == NULL) {
    bni_print_error("failed to read BAM header from %s", bam_path);
    sam_close(in);
    bni_index_destroy(&idx);
    free(default_index);
    return 1;
  }
  int status = 0;
  if (check_metadata(&idx, bam_path, hdr) != 0) {
    status = 1;
  }
  if (status == 0 && do_full) {
    BGZF *bgzf_fp = hts_get_bgzfp(in);
    if (bgzf_fp == NULL) {
      bni_print_error("failed to access BGZF handle");
      status = 1;
    } else if (full_check(&idx, in, hdr, bgzf_fp) != 0) {
      {
        status = 1;
      }
    }
  }
  if (status == 0) {
    printf("OK\n");
  }
  sam_hdr_destroy(hdr);
  sam_close(in);
  bni_index_destroy(&idx);
  free(default_index);
  return status;
}
