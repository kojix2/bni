#include "bni_internal.h"

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
#include <htslib/tbx.h> /* hts_get_bgzfp() is declared here in htslib */

typedef struct {
  bni_entry_t *data;
  size_t len;
  size_t cap;
} entry_vec_t;
typedef struct {
  char *data;
  size_t len;
  size_t cap;
} strbuf_t;
typedef struct {
  char *data;
  size_t cap;
} namebuf_t;

enum {
  INITIAL_ENTRY_CAPACITY = 4096,
  INITIAL_STRING_CAPACITY = 65536,
  BGZF_BLOCK_OFFSET_SHIFT = 16,
  FORMAT_BUFFER_SIZE = 64,
  OPT_NO_HEADER_CHECK = 1000,
};

static void usage_index(FILE *fp) {
  (void)fprintf(fp,
                "Usage:\n"
                "  bni index [options] <in.name.bam>\n\n"
                "Create <in.name.bam>.bni for a BAM sorted with samtools sort -N.\n"
                "BNIv2 stores one index entry per BGZF block that contains BAM record starts.\n\n"
                "Options:\n"
                "  -o, --output FILE          output index file [default: <in.bam>.bni]\n"
                "  -f, --force                overwrite an existing index\n"
                "  -@, --threads INT          decompression threads\n"
                "      --no-header-check      do not require @HD SO/SS tags\n"
                "  -v, --verbose              print progress summary\n"
                "  -h, --help                 show this help\n");
}

static void entry_vec_free(entry_vec_t *v) {
  free(v->data);
  v->data = NULL;
  v->len = v->cap = 0;
}
static void strbuf_free(strbuf_t *s) {
  free(s->data);
  s->data = NULL;
  s->len = s->cap = 0;
}
static void namebuf_free(namebuf_t *s) {
  free(s->data);
  s->data = NULL;
  s->cap = 0;
}

static int entry_vec_push(entry_vec_t *v, bni_entry_t e) {
  if (v->len == v->cap) {
    size_t new_cap = v->cap ? v->cap * 2 : INITIAL_ENTRY_CAPACITY;
    if (new_cap < v->cap || new_cap > SIZE_MAX / sizeof(bni_entry_t)) {
      bni_print_error("too many BGZF block entries for this platform");
      return -1;
    }
    bni_entry_t *p = (bni_entry_t *)realloc(v->data, new_cap * sizeof(bni_entry_t));
    if (p == NULL) {
      bni_print_error("out of memory while growing entry table");
      return -1;
    }
    v->data = p;
    v->cap = new_cap;
  }
  v->data[v->len++] = e;
  return 0;
}

static int strbuf_append_cstr(strbuf_t *s, const char *name, uint64_t *offset_out) {
  size_t n = strlen(name) + 1;
  if (s->len > UINT64_MAX) {
    return -1;
  }
  if (offset_out) {
    *offset_out = (uint64_t)s->len;
  }
  if (n > SIZE_MAX - s->len) {
    bni_print_error("string table is too large");
    return -1;
  }
  size_t need = s->len + n;
  if (need > s->cap) {
    size_t new_cap = s->cap ? s->cap * 2 : INITIAL_STRING_CAPACITY;
    while (new_cap < need) {
      if (new_cap > SIZE_MAX / 2) {
        new_cap = need;
        break;
      }
      new_cap *= 2;
    }
    char *p = (char *)realloc(s->data, new_cap);
    if (p == NULL) {
      bni_print_error("out of memory while growing string table");
      return -1;
    }
    s->data = p;
    s->cap = new_cap;
  }
  memcpy(s->data + s->len, name, n);
  s->len += n;
  return 0;
}

static int add_bgzf_block(entry_vec_t *entries, strbuf_t *strings, const char *first_name,
                          const char *last_name, uint64_t beg_voff, uint64_t end_voff,
                          uint32_t n_records) {
  uint64_t first_offset = 0;
  uint64_t last_offset = 0;
  if (first_name == NULL || last_name == NULL || n_records == 0) {
    bni_print_error("internal error while adding empty BGZF block entry");
    return -1;
  }
  if (strcmp(first_name, last_name) > 0) {
    bni_print_error("BAM is not lexicographically sorted inside a BGZF block near '%s' -> '%s'",
                    first_name, last_name);
    return -1;
  }
  if (beg_voff >= end_voff) {
    bni_print_error("internal error while adding empty virtual-offset range");
    return -1;
  }
  if (strbuf_append_cstr(strings, first_name, &first_offset) != 0) {
    return -1;
  }
  if (strbuf_append_cstr(strings, last_name, &last_offset) != 0) {
    return -1;
  }
  bni_entry_t e;
  memset(&e, 0, sizeof(e));
  e.first_name_offset = first_offset;
  e.last_name_offset = last_offset;
  e.beg_voff = beg_voff;
  e.end_voff = end_voff;
  e.n_records = n_records;
  return entry_vec_push(entries, e);
}

static int header_is_queryname_lex(sam_hdr_t *hdr, int no_header_check) {
  if (no_header_check) {
    return 0;
  }
  kstring_t so = {0, 0, NULL};
  kstring_t ss = {0, 0, NULL};
  int ok = 0;
  if (sam_hdr_find_tag_hd(hdr, "SO", &so) != 0 || so.s == NULL || strcmp(so.s, "queryname") != 0) {
    bni_print_error(
        "input header is not @HD SO:queryname; run: samtools sort -N -o out.name.bam in.bam");
    ok = -1;
    goto done;
  }
  if (sam_hdr_find_tag_hd(hdr, "SS", &ss) == 0 && ss.s != NULL) {
    if (strcmp(ss.s, "queryname:lexicographical") != 0) {
      bni_print_error("input header has @HD SS:%s, not SS:queryname:lexicographical", ss.s);
      ok = -1;
      goto done;
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
  if (text == NULL || len == SIZE_MAX) {
    return bni_fnv1a64("", 0);
  }
  return bni_fnv1a64(text, len);
}

static int namebuf_set(namebuf_t *dst, const char *src) {
  size_t n = strlen(src) + 1;
  if (n > dst->cap) {
    char *copy = (char *)realloc(dst->data, n);
    if (copy == NULL) {
      bni_print_error("out of memory while storing QNAME");
      return -1;
    }
    dst->data = copy;
    dst->cap = n;
  }
  memcpy(dst->data, src, n);
  return 0;
}

typedef struct {
  samFile *alignment;
  sam_hdr_t *header;
  BGZF *bgzf_fp;
} build_input_t;

typedef struct {
  entry_vec_t entries;
  strbuf_t strings;
  namebuf_t previous_name;
  namebuf_t block_first;
  uint64_t block_compressed_offset;
  uint64_t block_beg_voff;
  uint64_t total_records;
  uint64_t last_record_end_voff;
  uint32_t block_records;
  int have_block;
} index_scan_t;

static int resolve_build_output_path(const char *bam_path, const char *index_path, int force,
                                     char **default_out, const char **out_path) {
  *default_out = NULL;
  *out_path = index_path;
  if (*out_path == NULL) {
    *default_out = bni_default_index_path(bam_path);
    if (*default_out == NULL) {
      bni_print_error("out of memory while building output path");
      return -1;
    }
    *out_path = *default_out;
  }
  if (!force && bni_path_exists(*out_path)) {
    bni_print_error("%s already exists; use --force to overwrite", *out_path);
    return -1;
  }
  return 0;
}

static int read_bam_metadata(const char *bam_path, uint64_t *bam_size, int64_t *bam_mtime) {
  if (bni_file_metadata(bam_path, bam_size, bam_mtime) != 0) {
    bni_print_error("could not stat %s: %s", bam_path, strerror(errno));
    return -1;
  }
  return 0;
}

static int validate_build_input_format(samFile *alignment) {
  const htsFormat *fmt = hts_get_format(alignment);
  if (fmt != NULL && fmt->format != bam) {
    bni_print_error("input is not BAM; bni supports BGZF-compressed BAM only");
    return -1;
  }
  if (fmt != NULL && fmt->compression != bgzf) {
    bni_print_error("input is not BGZF-compressed BAM; bni cannot seek it safely");
    return -1;
  }
  return 0;
}

static int open_build_input(const char *bam_path, int threads, int no_header_check,
                            build_input_t *input) {
  input->alignment = sam_open(bam_path, "r");
  if (input->alignment == NULL) {
    bni_print_error("could not open %s", bam_path);
    return -1;
  }
  if (threads > 0 && hts_set_threads(input->alignment, threads) != 0) {
    bni_print_warning("failed to enable input threads");
  }
  input->header = sam_hdr_read(input->alignment);
  if (input->header == NULL) {
    bni_print_error("failed to read BAM header from %s", bam_path);
    return -1;
  }
  if (validate_build_input_format(input->alignment) != 0) {
    return -1;
  }
  if (header_is_queryname_lex(input->header, no_header_check) != 0) {
    return -1;
  }
  input->bgzf_fp = hts_get_bgzfp(input->alignment);
  if (input->bgzf_fp == NULL) {
    bni_print_error("failed to access BGZF handle");
    return -1;
  }
  return 0;
}

static void close_build_input(build_input_t *input) {
  sam_hdr_destroy(input->header);
  if (input->alignment != NULL) {
    sam_close(input->alignment);
  }
}

static void index_scan_destroy(index_scan_t *scan) {
  namebuf_free(&scan->previous_name);
  namebuf_free(&scan->block_first);
  entry_vec_free(&scan->entries);
  strbuf_free(&scan->strings);
}

static int check_qname_order(const index_scan_t *scan, const char *qname) {
  if (scan->previous_name.data == NULL || strcmp(scan->previous_name.data, qname) <= 0) {
    return 0;
  }
  bni_print_error(
      "BAM is not lexicographically queryname-sorted near '%s' -> '%s'; use samtools sort -N",
      scan->previous_name.data, qname);
  return -1;
}

static int start_scan_block(index_scan_t *scan, uint64_t record_beg_voff,
                            uint64_t record_compressed_offset, const char *qname) {
  scan->block_compressed_offset = record_compressed_offset;
  scan->block_beg_voff = record_beg_voff;
  if (namebuf_set(&scan->block_first, qname) != 0) {
    return -1;
  }
  scan->block_records = 1;
  scan->have_block = 1;
  return 0;
}

static int update_scan_block(index_scan_t *scan, uint64_t record_beg_voff, const char *qname) {
  uint64_t record_compressed_offset = record_beg_voff >> BGZF_BLOCK_OFFSET_SHIFT;
  if (!scan->have_block) {
    return start_scan_block(scan, record_beg_voff, record_compressed_offset, qname);
  }
  if (record_compressed_offset != scan->block_compressed_offset) {
    if (add_bgzf_block(&scan->entries, &scan->strings, scan->block_first.data,
                       scan->previous_name.data, scan->block_beg_voff, record_beg_voff,
                       scan->block_records) != 0) {
      return -1;
    }
    free(scan->previous_name.data); // NOLINT(clang-analyzer-unix.Malloc)
    scan->previous_name.data = NULL;
    scan->previous_name.cap = 0;
    return start_scan_block(scan, record_beg_voff, record_compressed_offset, qname);
  }
  if (scan->block_records == UINT32_MAX) {
    bni_print_error("too many records starting in one BGZF block");
    return -1;
  }
  scan->block_records++;
  return 0;
}

static int process_index_record(index_scan_t *scan, const char *qname, uint64_t record_beg_voff,
                                uint64_t record_end_voff) {
  if (qname == NULL || qname[0] == '\0') {
    bni_print_error("encountered record with empty QNAME");
    return -1;
  }
  if (check_qname_order(scan, qname) != 0) {
    return -1;
  }
  if (update_scan_block(scan, record_beg_voff, qname) != 0) {
    return -1;
  }
  if (namebuf_set(&scan->previous_name, qname) != 0) {
    return -1;
  }
  scan->last_record_end_voff = record_end_voff;
  scan->total_records++;
  return 0;
}

static int read_next_index_record(const build_input_t *input, bam1_t *record,
                                  uint64_t *record_end_voff) {
  int ret = sam_read1(input->alignment, input->header, record);
  if (ret < 0) {
    if (ret < -1) {
      bni_print_error("error while reading BAM records");
      return -1;
    }
    return 0;
  }
  int64_t tell_after = bgzf_tell(input->bgzf_fp);
  if (tell_after < 0) {
    bni_print_error("bgzf_tell failed after reading record");
    return -1;
  }
  *record_end_voff = (uint64_t)tell_after;
  return 1;
}

static int scan_bam_records(const build_input_t *input, index_scan_t *scan) {
  bam1_t *record = bam_init1();
  if (record == NULL) {
    bni_print_error("failed to allocate BAM record");
    return -1;
  }

  int status = -1;
  int64_t initial_voff = bgzf_tell(input->bgzf_fp);
  if (initial_voff < 0) {
    bni_print_error("bgzf_tell failed before reading first record");
    goto done;
  }
  uint64_t next_record_beg_voff = (uint64_t)initial_voff;
  for (;;) {
    uint64_t record_end_voff = 0;
    uint64_t record_beg_voff = next_record_beg_voff;
    int read_status = read_next_index_record(input, record, &record_end_voff);
    if (read_status < 0) {
      goto done;
    }
    if (read_status == 0) {
      break;
    }
    const char *qname = bam_get_qname(record);
    if (process_index_record(scan, qname, record_beg_voff, record_end_voff) != 0) {
      goto done;
    }
    next_record_beg_voff = record_end_voff;
  }
  if (scan->have_block && add_bgzf_block(&scan->entries, &scan->strings, scan->block_first.data,
                                         scan->previous_name.data, scan->block_beg_voff,
                                         scan->last_record_end_voff, scan->block_records) != 0) {
    goto done;
  }
  status = 0;

done:
  bam_destroy1(record);
  if (status != 0) {
    index_scan_destroy(scan);
  }
  return status;
}

static void init_index_header(bni_file_header_t *header, const index_scan_t *scan,
                              const build_input_t *input, uint64_t bam_size, int64_t bam_mtime) {
  memset(header, 0, sizeof(*header));
  header->version = BNI_FORMAT_VERSION;
  header->header_size = BNI_HEADER_SIZE;
  header->flags = BNI_FLAG_BGZF_BLOCKS;
  header->n_blocks = (uint64_t)scan->entries.len;
  header->n_records = scan->total_records;
  header->entries_offset = BNI_HEADER_SIZE;
  header->strings_offset =
      BNI_HEADER_SIZE + ((uint64_t)scan->entries.len * (uint64_t)BNI_ENTRY_SIZE);
  header->strings_size = (uint64_t)scan->strings.len;
  header->bam_size = bam_size;
  header->bam_mtime = bam_mtime;
  header->header_hash = header_hash64(input->header);
  header->sort_order = BNI_SORT_QUERYNAME_LEX;
  header->entry_size = BNI_ENTRY_SIZE;
}

static void copy_build_stats(const bni_file_header_t *header, bni_build_stats_t *stats) {
  if (stats != NULL) {
    stats->n_blocks = header->n_blocks;
    stats->n_records = header->n_records;
    stats->strings_size = header->strings_size;
  }
}

int bni_build_index(const char *bam_path, const char *index_path, const bni_build_options_t *opts,
                    bni_build_stats_t *stats) {
  int force = opts ? opts->force : 0;
  int threads = opts ? opts->threads : 0;
  int no_header_check = opts ? opts->no_header_check : 0;

  if (stats != NULL) {
    memset(stats, 0, sizeof(*stats));
  }
  if (bam_path == NULL || strcmp(bam_path, "-") == 0) {
    bni_print_error("indexing from stdin is not supported because BNI stores BGZF virtual offsets");
    return -1;
  }

  char *default_out = NULL;
  const char *out_path = NULL;
  uint64_t bam_size = 0;
  int64_t bam_mtime = 0;
  build_input_t input = {0};
  index_scan_t scan = {0};
  int status = -1;

  if (resolve_build_output_path(bam_path, index_path, force, &default_out, &out_path) != 0) {
    goto cleanup;
  }
  if (read_bam_metadata(bam_path, &bam_size, &bam_mtime) != 0) {
    goto cleanup;
  }
  if (open_build_input(bam_path, threads, no_header_check, &input) != 0) {
    goto cleanup;
  }
  if (scan_bam_records(&input, &scan) != 0) {
    goto cleanup;
  }

  bni_file_header_t header;
  init_index_header(&header, &scan, &input, bam_size, bam_mtime);
  const char *string_table = scan.strings.data ? scan.strings.data : "";
  if (bni_write_index_file(out_path, &header, scan.entries.data, string_table) != 0) {
    goto cleanup;
  }
  copy_build_stats(&header, stats);
  status = 0;

cleanup:
  index_scan_destroy(&scan);
  close_build_input(&input);
  free(default_out);
  return status;
}

int bni_cmd_index(int argc, char **argv) {
  const char *out_path_arg = NULL;
  int force = 0;
  int threads = 0;
  int no_header_check = 0;
  int verbose = 0;
  static const struct option long_opts[] = {
      {"output", required_argument, NULL, 'o'},
      {"force", no_argument, NULL, 'f'},
      {"threads", required_argument, NULL, '@'},
      {"no-header-check", no_argument, NULL, OPT_NO_HEADER_CHECK},
      {"verbose", no_argument, NULL, 'v'},
      {"help", no_argument, NULL, 'h'},
      {0, 0, 0, 0}};
  optind = 1;
  int c;
  while ((c = getopt_long(argc, argv, "o:f@::vh", long_opts, NULL)) != -1) {
    switch (c) {
    case 'o':
      out_path_arg = optarg;
      break;
    case 'f':
      force = 1;
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
    case 'v':
      verbose = 1;
      break;
    case OPT_NO_HEADER_CHECK:
      no_header_check = 1;
      break;
    case 'h':
      usage_index(stdout);
      return 0;
    default:
      usage_index(stderr);
      return 1;
    }
  }
  if (argc - optind != 1) {
    usage_index(stderr);
    return 1;
  }
  const char *bam_path = argv[optind];
  bni_build_options_t opts;
  memset(&opts, 0, sizeof(opts));
  opts.force = force;
  opts.threads = threads;
  opts.no_header_check = no_header_check;
  bni_build_stats_t stats;
  if (bni_build_index(bam_path, out_path_arg, &opts, &stats) != 0) {
    return 1;
  }
  if (verbose) {
    char bbuf[FORMAT_BUFFER_SIZE];
    char rbuf[FORMAT_BUFFER_SIZE];
    char sbuf[FORMAT_BUFFER_SIZE];
    char *default_out = NULL;
    const char *printed_out = out_path_arg;
    if (printed_out == NULL) {
      default_out = bni_default_index_path(bam_path);
      printed_out = default_out ? default_out : "(default index path)";
    }
    bni_format_u64(bbuf, sizeof(bbuf), stats.n_blocks);
    bni_format_u64(rbuf, sizeof(rbuf), stats.n_records);
    bni_format_u64(sbuf, sizeof(sbuf), stats.strings_size);
    (void)fprintf(stderr, "indexed %s: %s BGZF-block entries, %s records, %s string bytes -> %s\n",
                  bam_path, bbuf, rbuf, sbuf, printed_out);
    free(default_out);
  }
  return 0;
}
