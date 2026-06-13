#define _POSIX_C_SOURCE 200809L

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

static void usage_index(FILE *fp) {
  fprintf(fp, "Usage:\n"
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
    size_t new_cap = v->cap ? v->cap * 2 : 4096;
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
  if (s->len > UINT64_MAX)
    return -1;
  if (offset_out)
    *offset_out = (uint64_t)s->len;
  if (n > SIZE_MAX - s->len) {
    bni_print_error("string table is too large");
    return -1;
  }
  size_t need = s->len + n;
  if (need > s->cap) {
    size_t new_cap = s->cap ? s->cap * 2 : 65536;
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
  uint64_t first_offset = 0, last_offset = 0;
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
  if (strbuf_append_cstr(strings, first_name, &first_offset) != 0)
    return -1;
  if (strbuf_append_cstr(strings, last_name, &last_offset) != 0)
    return -1;
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
  if (no_header_check)
    return 0;
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
  if (text == NULL || len == SIZE_MAX)
    return bni_fnv1a64("", 0);
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

int bni_build_index(const char *bam_path, const char *index_path, const bni_build_options_t *opts,
                    bni_build_stats_t *stats) {
  int force = opts ? opts->force : 0;
  int threads = opts ? opts->threads : 0;
  int no_header_check = opts ? opts->no_header_check : 0;

  if (stats)
    memset(stats, 0, sizeof(*stats));
  if (bam_path == NULL || strcmp(bam_path, "-") == 0) {
    bni_print_error("indexing from stdin is not supported because BNI stores BGZF virtual offsets");
    return -1;
  }

  char *default_out = NULL;
  const char *out_path = index_path;
  if (out_path == NULL) {
    default_out = bni_default_index_path(bam_path);
    if (default_out == NULL) {
      bni_print_error("out of memory while building output path");
      return -1;
    }
    out_path = default_out;
  }
  if (!force && bni_path_exists(out_path)) {
    bni_print_error("%s already exists; use --force to overwrite", out_path);
    free(default_out);
    return -1;
  }
  uint64_t bam_size = 0;
  int64_t bam_mtime = 0;
  if (bni_file_metadata(bam_path, &bam_size, &bam_mtime) != 0) {
    bni_print_error("could not stat %s: %s", bam_path, strerror(errno));
    free(default_out);
    return -1;
  }

  samFile *fp = sam_open(bam_path, "r");
  if (fp == NULL) {
    bni_print_error("could not open %s", bam_path);
    free(default_out);
    return -1;
  }
  if (threads > 0 && hts_set_threads(fp, threads) != 0)
    bni_print_warning("failed to enable input threads");
  bam_hdr_t *hdr = sam_hdr_read(fp);
  if (hdr == NULL) {
    bni_print_error("failed to read BAM header from %s", bam_path);
    sam_close(fp);
    free(default_out);
    return -1;
  }
  const htsFormat *fmt = hts_get_format(fp);
  if (fmt != NULL && fmt->format != bam) {
    bni_print_error("input is not BAM; bni supports BGZF-compressed BAM only");
    sam_hdr_destroy(hdr);
    sam_close(fp);
    free(default_out);
    return -1;
  }
  if (fmt != NULL && fmt->compression != bgzf) {
    bni_print_error("input is not BGZF-compressed BAM; bni cannot seek it safely");
    sam_hdr_destroy(hdr);
    sam_close(fp);
    free(default_out);
    return -1;
  }
  if (header_is_queryname_lex(hdr, no_header_check) != 0) {
    sam_hdr_destroy(hdr);
    sam_close(fp);
    free(default_out);
    return -1;
  }
  BGZF *bgzf_fp = hts_get_bgzfp(fp);
  if (bgzf_fp == NULL) {
    bni_print_error("failed to access BGZF handle");
    sam_hdr_destroy(hdr);
    sam_close(fp);
    free(default_out);
    return -1;
  }

  entry_vec_t entries = {0, 0, 0};
  strbuf_t strings = {0, 0, 0};
  bam1_t *b = bam_init1();
  if (b == NULL) {
    bni_print_error("failed to allocate BAM record");
    sam_hdr_destroy(hdr);
    sam_close(fp);
    free(default_out);
    return -1;
  }

  namebuf_t prev_name = {0, 0};
  namebuf_t block_first = {0, 0};
  uint64_t block_coff = 0;
  uint64_t block_beg = 0;
  uint64_t total_records = 0;
  uint64_t last_rec_end = 0;
  uint32_t block_n = 0;
  int have_block = 0;
  int status = -1;

  for (;;) {
    int64_t tell_before = bgzf_tell(bgzf_fp);
    if (tell_before < 0) {
      bni_print_error("bgzf_tell failed before reading record");
      goto cleanup;
    }
    uint64_t rec_beg = (uint64_t)tell_before;
    int ret = sam_read1(fp, hdr, b);
    if (ret < 0) {
      if (ret < -1) {
        bni_print_error("error while reading BAM records");
        goto cleanup;
      }
      break;
    }
    int64_t tell_after = bgzf_tell(bgzf_fp);
    if (tell_after < 0) {
      bni_print_error("bgzf_tell failed after reading record");
      goto cleanup;
    }
    last_rec_end = (uint64_t)tell_after;
    const char *qname = bam_get_qname(b);
    if (qname == NULL || qname[0] == '\0') {
      bni_print_error("encountered record with empty QNAME");
      goto cleanup;
    }

    if (prev_name.data != NULL) {
      int cmp = strcmp(prev_name.data, qname);
      if (cmp > 0) {
        bni_print_error(
            "BAM is not lexicographically queryname-sorted near '%s' -> '%s'; use samtools sort -N",
            prev_name.data, qname);
        goto cleanup;
      }
    }

    uint64_t rec_coff = rec_beg >> 16;
    if (!have_block) {
      block_coff = rec_coff;
      block_beg = rec_beg;
      if (namebuf_set(&block_first, qname) != 0)
        goto cleanup;
      block_n = 1;
      have_block = 1;
    } else if (rec_coff != block_coff) {
      if (add_bgzf_block(&entries, &strings, block_first.data, prev_name.data, block_beg, rec_beg,
                         block_n) != 0)
        goto cleanup;
      block_coff = rec_coff;
      block_beg = rec_beg;
      if (namebuf_set(&block_first, qname) != 0)
        goto cleanup;
      block_n = 1;
    } else {
      if (block_n == UINT32_MAX) {
        bni_print_error("too many records starting in one BGZF block");
        goto cleanup;
      }
      block_n++;
    }

    if (namebuf_set(&prev_name, qname) != 0)
      goto cleanup;
    total_records++;
  }

  if (have_block) {
    if (add_bgzf_block(&entries, &strings, block_first.data, prev_name.data, block_beg,
                       last_rec_end, block_n) != 0)
      goto cleanup;
  }

  bni_file_header_t header;
  memset(&header, 0, sizeof(header));
  header.version = BNI_FORMAT_VERSION;
  header.header_size = BNI_HEADER_SIZE;
  header.flags = BNI_FLAG_BGZF_BLOCKS;
  header.n_blocks = (uint64_t)entries.len;
  header.n_records = total_records;
  header.entries_offset = BNI_HEADER_SIZE;
  header.strings_offset = BNI_HEADER_SIZE + (uint64_t)entries.len * (uint64_t)BNI_ENTRY_SIZE;
  header.strings_size = (uint64_t)strings.len;
  header.bam_size = bam_size;
  header.bam_mtime = bam_mtime;
  header.header_hash = header_hash64(hdr);
  header.sort_order = BNI_SORT_QUERYNAME_LEX;
  header.entry_size = BNI_ENTRY_SIZE;
  if (bni_write_index_file(out_path, &header, entries.data, strings.data ? strings.data : "") != 0)
    goto cleanup;
  if (stats) {
    stats->n_blocks = header.n_blocks;
    stats->n_records = header.n_records;
    stats->strings_size = header.strings_size;
  }
  status = 0;

cleanup:
  namebuf_free(&prev_name);
  namebuf_free(&block_first);
  bam_destroy1(b);
  entry_vec_free(&entries);
  strbuf_free(&strings);
  sam_hdr_destroy(hdr);
  sam_close(fp);
  free(default_out);
  return status;
}

int bni_cmd_index(int argc, char **argv) {
  const char *out_path_arg = NULL;
  int force = 0, threads = 0, no_header_check = 0, verbose = 0;
  static const struct option long_opts[] = {{"output", required_argument, NULL, 'o'},
                                            {"force", no_argument, NULL, 'f'},
                                            {"threads", required_argument, NULL, '@'},
                                            {"no-header-check", no_argument, NULL, 1000},
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
    case 1000:
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
  if (bni_build_index(bam_path, out_path_arg, &opts, &stats) != 0)
    return 1;
  if (verbose) {
    char bbuf[64], rbuf[64], sbuf[64];
    char *default_out = NULL;
    const char *printed_out = out_path_arg;
    if (printed_out == NULL) {
      default_out = bni_default_index_path(bam_path);
      printed_out = default_out ? default_out : "(default index path)";
    }
    bni_format_u64(bbuf, sizeof(bbuf), stats.n_blocks);
    bni_format_u64(rbuf, sizeof(rbuf), stats.n_records);
    bni_format_u64(sbuf, sizeof(sbuf), stats.strings_size);
    fprintf(stderr, "indexed %s: %s BGZF-block entries, %s records, %s string bytes -> %s\n",
            bam_path, bbuf, rbuf, sbuf, printed_out);
    free(default_out);
  }
  return 0;
}
