#include "bni_internal.h"

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
#include <htslib/tbx.h> /* hts_get_bgzfp() is declared here in htslib */

typedef enum { BNI_OUT_BAM, BNI_OUT_SAM, BNI_OUT_CRAM } out_format_t;

enum {
  OPT_NO_HEADER = 1000,
  OPT_WITH_HEADER = 1001,
  OPT_MISSING_OK = 1002,
  OPT_LIST_MISSING = 1003,
  OPT_IGNORE_METADATA = 1004,
  INITIAL_NAME_REQUEST_CAPACITY = 1024,
  FORMAT_BUFFER_SIZE = 64,
};

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
  if (strcmp(s, "bam") == 0) {
    *fmt = BNI_OUT_BAM;
    return 0;
  }
  if (strcmp(s, "sam") == 0) {
    *fmt = BNI_OUT_SAM;
    return 0;
  }
  if (strcmp(s, "cram") == 0) {
    *fmt = BNI_OUT_CRAM;
    return 0;
  }
  return -1;
}

static out_format_t infer_output_format(const char *path) {
  if (path != NULL) {
    if (bni_has_suffix(path, ".sam")) {
      return BNI_OUT_SAM;
    }
    if (bni_has_suffix(path, ".cram")) {
      return BNI_OUT_CRAM;
    }
  }
  return BNI_OUT_BAM;
}

static const char *mode_for_format(out_format_t fmt) {
  switch (fmt) {
  case BNI_OUT_BAM:
    return "wb";
  case BNI_OUT_SAM:
    return "w";
  case BNI_OUT_CRAM:
    return "wc";
  }
  return "wb";
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
  if (bni_file_metadata(bam_path, &bam_size, &bam_mtime) != 0) {
    bni_print_error("could not stat %s: %s", bam_path, strerror(errno));
    return -1;
  }
  if (idx->header.bam_size != bam_size) {
    bni_print_error(
        "BAM size differs from index metadata; rebuild the BNI index or use --ignore-metadata");
    return -1;
  }
  if (idx->header.bam_mtime != bam_mtime) {
    bni_print_error(
        "BAM mtime differs from index metadata; rebuild the BNI index or use --ignore-metadata");
    return -1;
  }
  uint64_t hh = header_hash64(hdr);
  if (idx->header.header_hash != hh) {
    bni_print_error("BAM header hash differs from index metadata; rebuild the BNI index or use "
                    "--ignore-metadata");
    return -1;
  }
  return 0;
}

static char *trim_line(char *line) {
  if (line == NULL) {
    return NULL;
  }
  size_t n = strlen(line);
  while (n > 0 && (line[n - 1] == '\n' || line[n - 1] == '\r')) {
    line[--n] = '\0';
  }
  return line;
}

struct bni_reader_t {
  bni_index_t idx;
  samFile *in;
  sam_hdr_t *hdr;
  BGZF *bgzf_fp;
  bam1_t *record;
};

int bni_reader_open(const char *bam_path, const char *index_path, const bni_reader_options_t *opts,
                    bni_reader_t **reader_out) {
  int threads = opts ? opts->threads : 0;
  int ignore_metadata = opts ? opts->ignore_metadata : 0;

  if (reader_out == NULL) {
    return -1;
  }
  *reader_out = NULL;
  if (bam_path == NULL || strcmp(bam_path, "-") == 0) {
    bni_print_error("random access from stdin is not supported");
    return -1;
  }

  char *default_index = NULL;
  const char *resolved_index_path = index_path;
  if (resolved_index_path == NULL) {
    default_index = bni_default_index_path(bam_path);
    if (default_index == NULL) {
      bni_print_error("out of memory while building default index path");
      return -1;
    }
    resolved_index_path = default_index;
  }

  bni_reader_t *reader = (bni_reader_t *)calloc(1, sizeof(*reader));
  if (reader == NULL) {
    bni_print_error("out of memory while allocating reader");
    free(default_index);
    return -1;
  }
  if (bni_load_index_file(resolved_index_path, &reader->idx) != 0) {
    goto fail;
  }
  reader->in = sam_open(bam_path, "r");
  if (reader->in == NULL) {
    bni_print_error("could not open %s", bam_path);
    goto fail;
  }
  if (threads > 0 && hts_set_threads(reader->in, threads) != 0) {
    bni_print_warning("failed to enable input threads");
  }
  reader->hdr = sam_hdr_read(reader->in);
  if (reader->hdr == NULL) {
    bni_print_error("failed to read BAM header from %s", bam_path);
    goto fail;
  }
  const htsFormat *fmt = hts_get_format(reader->in);
  if (fmt != NULL && fmt->format != bam) {
    bni_print_error("input is not BAM; bni supports BGZF-compressed BAM only");
    goto fail;
  }
  if (!ignore_metadata && check_metadata(&reader->idx, bam_path, reader->hdr) != 0) {
    goto fail;
  }
  reader->bgzf_fp = hts_get_bgzfp(reader->in);
  if (reader->bgzf_fp == NULL) {
    bni_print_error("failed to access BGZF handle");
    goto fail;
  }
  reader->record = bam_init1();
  if (reader->record == NULL) {
    bni_print_error("failed to allocate BAM record");
    goto fail;
  }

  free(default_index);
  *reader_out = reader;
  return 0;

fail:
  free(default_index);
  bni_reader_close(reader);
  return -1;
}

void bni_reader_close(bni_reader_t *reader) {
  if (reader == NULL) {
    return;
  }
  bam_destroy1(reader->record);
  sam_hdr_destroy(reader->hdr);
  if (reader->in) {
    sam_close(reader->in);
  }
  bni_index_destroy(&reader->idx);
  free(reader);
}

const bni_index_t *bni_reader_index(const bni_reader_t *reader) {
  return reader ? &reader->idx : NULL;
}

const sam_hdr_t *bni_reader_header(const bni_reader_t *reader) {
  return reader ? reader->hdr : NULL;
}

int bni_reader_fetch(bni_reader_t *reader, const char *name, bni_record_callback callback,
                     void *user, uint32_t *n_records_out) {
  if (n_records_out) {
    *n_records_out = 0;
  }
  if (reader == NULL || name == NULL) {
    return -1;
  }

  const bni_entry_t *entry = bni_find_entry(&reader->idx, name);
  if (entry == NULL) {
    return 1;
  }

  if (bgzf_seek(reader->bgzf_fp, (int64_t)entry->beg_voff, SEEK_SET) < 0) {
    bni_print_error("bgzf_seek failed for QNAME '%s'", name);
    return -1;
  }

  uint32_t seen = 0;
  for (;;) {
    int ret = sam_read1(reader->in, reader->hdr, reader->record);
    if (ret < 0) {
      if (ret < -1) {
        bni_print_error("error while reading BAM while fetching QNAME '%s'", name);
        return -1;
      }
      break;
    }
    const char *qname = bam_get_qname(reader->record);
    if (qname == NULL) {
      bni_print_error("encountered record with NULL QNAME while fetching '%s'", name);
      return -1;
    }
    int cmp = strcmp(qname, name);
    if (cmp < 0) {
      continue;
    }
    if (cmp > 0) {
      break;
    }
    if (callback != NULL && callback(reader->record, reader->hdr, user) != 0) {
      return -1;
    }
    if (seen == UINT32_MAX) {
      bni_print_error("too many records matched QNAME '%s'", name);
      return -1;
    }
    seen++;
  }

  if (seen == 0) {
    return 1;
  }
  if (n_records_out) {
    *n_records_out = seen;
  }
  return 0;
}

typedef struct {
  samFile *out;
} write_context_t;

typedef struct {
  char *name;
  uint64_t entry_index;
  int found;
} name_request_t;

typedef struct {
  name_request_t *data;
  size_t len;
  size_t cap;
} name_request_vec_t;

static int write_record_callback(const bam1_t *record, const sam_hdr_t *hdr, void *user) {
  write_context_t *ctx = (write_context_t *)user;
  if (sam_write1(ctx->out, hdr, record) < 0) {
    bni_print_error("failed writing output record");
    return -1;
  }
  return 0;
}

static void name_request_vec_destroy(name_request_vec_t *v) {
  if (v == NULL) {
    return;
  }
  for (size_t i = 0; i < v->len; ++i) {
    free(v->data[i].name);
  }
  free(v->data);
  v->data = NULL;
  v->len = v->cap = 0;
}

static int name_request_vec_push(name_request_vec_t *v, char *name, uint64_t entry_index) {
  if (v->len == v->cap) {
    size_t new_cap = v->cap ? v->cap * 2 : INITIAL_NAME_REQUEST_CAPACITY;
    if (new_cap < v->cap || new_cap > SIZE_MAX / sizeof(name_request_t)) {
      bni_print_error("too many requested read names for this platform");
      return -1;
    }
    name_request_t *p = (name_request_t *)realloc(v->data, new_cap * sizeof(*v->data));
    if (p == NULL) {
      bni_print_error("out of memory while loading requested read names");
      return -1;
    }
    v->data = p;
    v->cap = new_cap;
  }
  v->data[v->len].name = name;
  v->data[v->len].entry_index = entry_index;
  v->data[v->len].found = 0;
  v->len++;
  return 0;
}

static int compare_name_requests(const void *a, const void *b) {
  const name_request_t *ra = (const name_request_t *)a;
  const name_request_t *rb = (const name_request_t *)b;
  if (ra->entry_index < rb->entry_index) {
    return -1;
  }
  if (ra->entry_index > rb->entry_index) {
    return 1;
  }
  return strcmp(ra->name, rb->name);
}

static void dedupe_name_requests(name_request_vec_t *v) {
  if (v->len == 0) {
    return;
  }
  size_t out = 1;
  for (size_t i = 1; i < v->len; ++i) {
    name_request_t *prev = &v->data[out - 1];
    name_request_t *cur = &v->data[i];
    if (prev->entry_index == cur->entry_index && strcmp(prev->name, cur->name) == 0) {
      free(cur->name);
      continue;
    }
    if (out != i) {
      v->data[out] = *cur;
    }
    out++;
  }
  v->len = out;
}

static int load_name_requests(FILE *nf, bni_reader_t *reader, name_request_vec_t *requests,
                              int list_missing, uint64_t *missing_out) {
  char *line = NULL;
  size_t cap = 0;
  ssize_t nread;
  int status = 0;
  while ((nread = getline(&line, &cap, nf)) >= 0) {
    (void)nread;
    char *trimmed = trim_line(line);
    if (trimmed[0] == '\0') {
      continue;
    }

    const bni_entry_t *entry = bni_find_entry(&reader->idx, trimmed);
    if (entry == NULL) {
      if (list_missing) {
        fprintf(stderr, "%s\n", trimmed);
      }
      (*missing_out)++;
      continue;
    }

    char *name = bni_strdup(trimmed);
    if (name == NULL) {
      bni_print_error("out of memory while loading requested read names");
      status = -1;
      break;
    }
    uint64_t entry_index = (uint64_t)(entry - reader->idx.entries);
    if (name_request_vec_push(requests, name, entry_index) != 0) {
      free(name);
      status = -1;
      break;
    }
  }
  free(line);
  return status;
}

static int fetch_name_requests(bni_reader_t *reader, name_request_vec_t *requests,
                               write_context_t *write_ctx, int list_missing,
                               uint64_t *missing_out) {
  if (requests->len == 0) {
    return 0;
  }
  qsort(requests->data, requests->len, sizeof(*requests->data), compare_name_requests);
  dedupe_name_requests(requests);

  size_t group_beg = 0;
  while (group_beg < requests->len) {
    uint64_t entry_index = requests->data[group_beg].entry_index;
    size_t group_end = group_beg + 1;
    while (group_end < requests->len && requests->data[group_end].entry_index == entry_index) {
      group_end++;
    }

    const bni_entry_t *entry = &reader->idx.entries[entry_index];
    if (bgzf_seek(reader->bgzf_fp, (int64_t)entry->beg_voff, SEEK_SET) < 0) {
      bni_print_error("bgzf_seek failed for indexed name batch");
      return -1;
    }

    size_t target = group_beg;
    const char *last_target = requests->data[group_end - 1].name;
    for (;;) {
      int ret = sam_read1(reader->in, reader->hdr, reader->record);
      if (ret < 0) {
        if (ret < -1) {
          bni_print_error("error while reading BAM while fetching requested names");
          return -1;
        }
        break;
      }
      const char *qname = bam_get_qname(reader->record);
      if (qname == NULL) {
        bni_print_error("encountered record with NULL QNAME while fetching requested names");
        return -1;
      }

      if (strcmp(qname, last_target) > 0) {
        break;
      }
      while (target < group_end && strcmp(requests->data[target].name, qname) < 0) {
        if (!requests->data[target].found) {
          if (list_missing) {
            fprintf(stderr, "%s\n", requests->data[target].name);
          }
          (*missing_out)++;
        }
        target++;
      }
      if (target == group_end) {
        break;
      }
      if (strcmp(qname, requests->data[target].name) == 0) {
        requests->data[target].found = 1;
        if (write_record_callback(reader->record, reader->hdr, write_ctx) != 0) {
          return -1;
        }
      }
    }

    while (target < group_end) {
      if (!requests->data[target].found) {
        if (list_missing) {
          fprintf(stderr, "%s\n", requests->data[target].name);
        }
        (*missing_out)++;
      }
      target++;
    }
    group_beg = group_end;
  }
  return 0;
}

int bni_cmd_get(int argc, char **argv) {
  const char *index_path_arg = NULL;
  const char *out_path = NULL;
  const char *names_path = NULL;
  const char *fmt_arg = NULL;
  int write_header = 1;
  int missing_ok = 0;
  int list_missing = 0;
  int ignore_metadata = 0;
  int threads = 0;
  static const struct option long_opts[] = {
      {"index", required_argument, NULL, 'i'},
      {"output", required_argument, NULL, 'o'},
      {"output-format", required_argument, NULL, 'O'},
      {"name-file", required_argument, NULL, 'f'},
      {"threads", required_argument, NULL, '@'},
      {"no-header", no_argument, NULL, OPT_NO_HEADER},
      {"with-header", no_argument, NULL, OPT_WITH_HEADER},
      {"missing-ok", no_argument, NULL, OPT_MISSING_OK},
      {"list-missing", no_argument, NULL, OPT_LIST_MISSING},
      {"ignore-metadata", no_argument, NULL, OPT_IGNORE_METADATA},
      {"help", no_argument, NULL, 'h'},
      {0, 0, 0, 0}};
  optind = 1;
  int c;
  while ((c = getopt_long(argc, argv, "i:o:O:f:@::h", long_opts, NULL)) != -1) {
    switch (c) {
    case 'i':
      index_path_arg = optarg;
      break;
    case 'o':
      out_path = optarg;
      break;
    case 'O':
      fmt_arg = optarg;
      break;
    case 'f':
      names_path = optarg;
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
    case OPT_NO_HEADER:
      write_header = 0;
      break;
    case OPT_WITH_HEADER:
      write_header = 1;
      break;
    case OPT_MISSING_OK:
      missing_ok = 1;
      break;
    case OPT_LIST_MISSING:
      list_missing = 1;
      break;
    case OPT_IGNORE_METADATA:
      ignore_metadata = 1;
      break;
    case 'h':
      usage_get(stdout);
      return 0;
    default:
      usage_get(stderr);
      return 1;
    }
  }
  const char *bam_path = NULL;
  const char *single_name = NULL;
  if (names_path != NULL) {
    if (argc - optind != 1) {
      usage_get(stderr);
      return 1;
    }
    bam_path = argv[optind];
  } else {
    if (argc - optind != 2) {
      usage_get(stderr);
      return 1;
    }
    bam_path = argv[optind];
    single_name = argv[optind + 1];
  }
  if (strcmp(bam_path, "-") == 0) {
    bni_print_error("random access from stdin is not supported");
    return 1;
  }
  out_format_t out_fmt = infer_output_format(out_path);
  if (fmt_arg != NULL && parse_output_format(fmt_arg, &out_fmt) != 0) {
    bni_print_error("unknown output format '%s'; expected bam, sam, or cram", fmt_arg);
    return 1;
  }
  if (!write_header && out_fmt != BNI_OUT_SAM) {
    bni_print_error("--no-header is only supported for SAM output");
    return 1;
  }
  bni_reader_options_t reader_opts;
  memset(&reader_opts, 0, sizeof(reader_opts));
  reader_opts.threads = threads;
  reader_opts.ignore_metadata = ignore_metadata;
  bni_reader_t *reader = NULL;
  if (bni_reader_open(bam_path, index_path_arg, &reader_opts, &reader) != 0) {
    return 1;
  }
  samFile *out = sam_open(out_path ? out_path : "-", mode_for_format(out_fmt));
  if (out == NULL) {
    bni_print_error("could not open output");
    bni_reader_close(reader);
    return 1;
  }
  if (threads > 0 && hts_set_threads(out, threads) != 0) {
    bni_print_warning("failed to enable output threads");
  }
  if (write_header && sam_hdr_write(out, bni_reader_header(reader)) != 0) {
    bni_print_error("failed to write output header");
    sam_close(out);
    bni_reader_close(reader);
    return 1;
  }
  write_context_t write_ctx;
  write_ctx.out = out;
  int status = 0;
  uint64_t missing = 0;
  if (names_path != NULL) {
    int close_nf = strcmp(names_path, "-") != 0;
    FILE *nf = close_nf ? fopen(names_path, "r") : stdin;
    if (nf == NULL) {
      bni_print_error("could not open name file %s: %s", names_path, strerror(errno));
      status = 1;
    } else {
      name_request_vec_t requests = {0, 0, 0};
      if (load_name_requests(nf, reader, &requests, list_missing, &missing) != 0 ||
          fetch_name_requests(reader, &requests, &write_ctx, list_missing, &missing) != 0) {
        status = 1;
      }
      name_request_vec_destroy(&requests);
      if (close_nf) {
        fclose(nf);
      }
    }
  } else {
    int ret = bni_reader_fetch(reader, single_name, write_record_callback, &write_ctx, NULL);
    if (ret > 0) {
      if (list_missing) {
        fprintf(stderr, "%s\n", single_name);
      }
      missing++;
    } else if (ret < 0) {
      {
        status = 1;
      }
    }
  }
  if (sam_close(out) != 0) {
    bni_print_error("failed closing output");
    status = 1;
  }
  bni_reader_close(reader);
  if (status == 0 && missing > 0 && !missing_ok) {
    char mbuf[FORMAT_BUFFER_SIZE];
    bni_format_u64(mbuf, sizeof(mbuf), missing);
    bni_print_error("%s requested read name(s) were not found; use --missing-ok to ignore", mbuf);
    status = 1;
  }
  return status;
}
