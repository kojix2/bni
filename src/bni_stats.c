#include "bni_internal.h"

#include <getopt.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

enum {
  TIME_BUFFER_SIZE = 64,
  FORMAT_BUFFER_SIZE = 64,
};

static void usage_stats(FILE *fp) {
  fprintf(fp, "Usage:\n"
              "  bni stats [options] <in.name.bam>\n\n"
              "Options:\n"
              "  -i, --index FILE       index file [default: <in.bam>.bni]\n"
              "  -h, --help             show this help\n");
}

static void print_mtime(int64_t t) {
  time_t tt = (time_t)t;
  struct tm tm_buf;
  struct tm *tm = localtime_r(&tt, &tm_buf);
  if (tm == NULL) {
    printf("bam_mtime: %" PRId64 "\n", t);
    return;
  }
  char buf[TIME_BUFFER_SIZE];
  if (strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S%z", tm) == 0) {
    printf("bam_mtime: %" PRId64 "\n", t);
    return;
  }
  printf("bam_mtime: %s\n", buf);
}

int bni_cmd_stats(int argc, char **argv) {
  const char *index_path_arg = NULL;
  static const struct option long_opts[] = {
      {"index", required_argument, NULL, 'i'}, {"help", no_argument, NULL, 'h'}, {0, 0, 0, 0}};
  optind = 1;
  int c;
  while ((c = getopt_long(argc, argv, "i:h", long_opts, NULL)) != -1) {
    switch (c) {
    case 'i':
      index_path_arg = optarg;
      break;
    case 'h':
      usage_stats(stdout);
      return 0;
    default:
      usage_stats(stderr);
      return 1;
    }
  }
  if (argc - optind != 1) {
    usage_stats(stderr);
    return 1;
  }
  const char *bam_path = argv[optind];
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
  char block_buf[FORMAT_BUFFER_SIZE];
  char rec_buf[FORMAT_BUFFER_SIZE];
  char str_buf[FORMAT_BUFFER_SIZE];
  char bam_buf[FORMAT_BUFFER_SIZE];
  bni_format_u64(block_buf, sizeof(block_buf), idx.header.n_blocks);
  bni_format_u64(rec_buf, sizeof(rec_buf), idx.header.n_records);
  bni_format_u64(str_buf, sizeof(str_buf), idx.header.strings_size);
  bni_format_u64(bam_buf, sizeof(bam_buf), idx.header.bam_size);
  printf("index: %s\n", index_path);
  printf("bam: %s\n", bam_path);
  printf("format: BNIv%u\n", idx.header.version);
  printf("mode: bgzf-block-name-ranges\n");
  printf("sort_order: queryname:lexicographical\n");
  printf("blocks: %s\n", block_buf);
  printf("records: %s\n", rec_buf);
  if (idx.header.n_blocks > 0) {
    printf("avg_records_per_block: %.6f\n",
           (double)idx.header.n_records / (double)idx.header.n_blocks);
  } else {
    printf("avg_records_per_block: 0\n");
  }
  printf("entry_size: %u\n", idx.header.entry_size);
  printf("string_table_bytes: %s\n", str_buf);
  printf("bam_size: %s\n", bam_buf);
  print_mtime(idx.header.bam_mtime);
  printf("header_hash_fnv1a64: %016" PRIx64 "\n", idx.header.header_hash);
  bni_index_destroy(&idx);
  free(default_index);
  return 0;
}
