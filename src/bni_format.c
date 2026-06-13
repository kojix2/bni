#include "bni_internal.h"

#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#define BNI_ENTRY_WRITE_CHUNK 8192u

enum {
  BNI_BITS_PER_BYTE = 8,
  BNI_U64_BYTES = 8,
  BNI_U8_MASK = 0xffU,

  BNI_HEADER_MAGIC_OFFSET = 0,
  BNI_HEADER_VERSION_OFFSET = 4,
  BNI_HEADER_SIZE_OFFSET = 8,
  BNI_HEADER_FLAGS_OFFSET = 12,
  BNI_HEADER_N_BLOCKS_OFFSET = 16,
  BNI_HEADER_N_RECORDS_OFFSET = 24,
  BNI_HEADER_ENTRIES_OFFSET_OFFSET = 32,
  BNI_HEADER_STRINGS_OFFSET_OFFSET = 40,
  BNI_HEADER_STRINGS_SIZE_OFFSET = 48,
  BNI_HEADER_BAM_SIZE_OFFSET = 56,
  BNI_HEADER_BAM_MTIME_OFFSET = 64,
  BNI_HEADER_HASH_OFFSET = 72,
  BNI_HEADER_SORT_ORDER_OFFSET = 80,
  BNI_HEADER_ENTRY_SIZE_OFFSET = 84,

  BNI_ENTRY_FIRST_NAME_OFFSET = 0,
  BNI_ENTRY_LAST_NAME_OFFSET = 8,
  BNI_ENTRY_BEG_VOFF_OFFSET = 16,
  BNI_ENTRY_END_VOFF_OFFSET = 24,
  BNI_ENTRY_N_RECORDS_OFFSET = 32,
  BNI_ENTRY_RESERVED_OFFSET = 36,
};

static void put_u32le(unsigned char *p, uint32_t v) {
  p[0] = (unsigned char)(v & BNI_U8_MASK);
  p[1] = (unsigned char)((v >> BNI_BITS_PER_BYTE) & BNI_U8_MASK);
  p[2] = (unsigned char)((v >> (2 * BNI_BITS_PER_BYTE)) & BNI_U8_MASK);
  p[3] = (unsigned char)((v >> (3 * BNI_BITS_PER_BYTE)) & BNI_U8_MASK);
}

static void put_u64le(unsigned char *p, uint64_t v) {
  for (int i = 0; i < BNI_U64_BYTES; ++i) {
    p[i] = (unsigned char)((v >> (BNI_BITS_PER_BYTE * i)) & BNI_U8_MASK);
  }
}

static uint32_t get_u32le(const unsigned char *p) {
  return ((uint32_t)p[0]) | ((uint32_t)p[1] << BNI_BITS_PER_BYTE) |
         ((uint32_t)p[2] << (2 * BNI_BITS_PER_BYTE)) | ((uint32_t)p[3] << (3 * BNI_BITS_PER_BYTE));
}

static uint64_t get_u64le(const unsigned char *p) {
  uint64_t v = 0;
  for (int i = BNI_U64_BYTES - 1; i >= 0; --i) {
    v <<= BNI_BITS_PER_BYTE;
    v |= (uint64_t)p[i];
  }
  return v;
}

static int write_exact(FILE *fp, const void *buf, size_t n) {
  return fwrite(buf, 1, n, fp) == n ? 0 : -1;
}

static int read_exact(FILE *fp, void *buf, size_t n) { return fread(buf, 1, n, fp) == n ? 0 : -1; }

static void close_ignoring_error(FILE *fp) { (void)fclose(fp); }

static void encode_header(unsigned char out[BNI_HEADER_SIZE], const bni_file_header_t *h) {
  memset(out, 0, BNI_HEADER_SIZE);
  out[BNI_HEADER_MAGIC_OFFSET] = BNI_MAGIC0;
  out[BNI_HEADER_MAGIC_OFFSET + 1] = BNI_MAGIC1;
  out[BNI_HEADER_MAGIC_OFFSET + 2] = BNI_MAGIC2;
  out[BNI_HEADER_MAGIC_OFFSET + 3] = BNI_MAGIC3;
  put_u32le(out + BNI_HEADER_VERSION_OFFSET, h->version);
  put_u32le(out + BNI_HEADER_SIZE_OFFSET, h->header_size);
  put_u32le(out + BNI_HEADER_FLAGS_OFFSET, h->flags);
  put_u64le(out + BNI_HEADER_N_BLOCKS_OFFSET, h->n_blocks);
  put_u64le(out + BNI_HEADER_N_RECORDS_OFFSET, h->n_records);
  put_u64le(out + BNI_HEADER_ENTRIES_OFFSET_OFFSET, h->entries_offset);
  put_u64le(out + BNI_HEADER_STRINGS_OFFSET_OFFSET, h->strings_offset);
  put_u64le(out + BNI_HEADER_STRINGS_SIZE_OFFSET, h->strings_size);
  put_u64le(out + BNI_HEADER_BAM_SIZE_OFFSET, h->bam_size);
  put_u64le(out + BNI_HEADER_BAM_MTIME_OFFSET, (uint64_t)h->bam_mtime);
  put_u64le(out + BNI_HEADER_HASH_OFFSET, h->header_hash);
  put_u32le(out + BNI_HEADER_SORT_ORDER_OFFSET, h->sort_order);
  put_u32le(out + BNI_HEADER_ENTRY_SIZE_OFFSET, h->entry_size);
}

static int decode_header(const unsigned char in[BNI_HEADER_SIZE], bni_file_header_t *h) {
  if (in[BNI_HEADER_MAGIC_OFFSET] != BNI_MAGIC0 || in[BNI_HEADER_MAGIC_OFFSET + 1] != BNI_MAGIC1 ||
      in[BNI_HEADER_MAGIC_OFFSET + 2] != BNI_MAGIC2 ||
      in[BNI_HEADER_MAGIC_OFFSET + 3] != BNI_MAGIC3) {
    return -1;
  }
  memset(h, 0, sizeof(*h));
  h->version = get_u32le(in + BNI_HEADER_VERSION_OFFSET);
  h->header_size = get_u32le(in + BNI_HEADER_SIZE_OFFSET);
  h->flags = get_u32le(in + BNI_HEADER_FLAGS_OFFSET);
  h->n_blocks = get_u64le(in + BNI_HEADER_N_BLOCKS_OFFSET);
  h->n_records = get_u64le(in + BNI_HEADER_N_RECORDS_OFFSET);
  h->entries_offset = get_u64le(in + BNI_HEADER_ENTRIES_OFFSET_OFFSET);
  h->strings_offset = get_u64le(in + BNI_HEADER_STRINGS_OFFSET_OFFSET);
  h->strings_size = get_u64le(in + BNI_HEADER_STRINGS_SIZE_OFFSET);
  h->bam_size = get_u64le(in + BNI_HEADER_BAM_SIZE_OFFSET);
  h->bam_mtime = (int64_t)get_u64le(in + BNI_HEADER_BAM_MTIME_OFFSET);
  h->header_hash = get_u64le(in + BNI_HEADER_HASH_OFFSET);
  h->sort_order = get_u32le(in + BNI_HEADER_SORT_ORDER_OFFSET);
  h->entry_size = get_u32le(in + BNI_HEADER_ENTRY_SIZE_OFFSET);
  return 0;
}

static void encode_entry(unsigned char out[BNI_ENTRY_SIZE], const bni_entry_t *e) {
  memset(out, 0, BNI_ENTRY_SIZE);
  put_u64le(out + BNI_ENTRY_FIRST_NAME_OFFSET, e->first_name_offset);
  put_u64le(out + BNI_ENTRY_LAST_NAME_OFFSET, e->last_name_offset);
  put_u64le(out + BNI_ENTRY_BEG_VOFF_OFFSET, e->beg_voff);
  put_u64le(out + BNI_ENTRY_END_VOFF_OFFSET, e->end_voff);
  put_u32le(out + BNI_ENTRY_N_RECORDS_OFFSET, e->n_records);
  put_u32le(out + BNI_ENTRY_RESERVED_OFFSET, e->reserved);
}

static void decode_entry(const unsigned char in[BNI_ENTRY_SIZE], bni_entry_t *e) {
  e->first_name_offset = get_u64le(in + BNI_ENTRY_FIRST_NAME_OFFSET);
  e->last_name_offset = get_u64le(in + BNI_ENTRY_LAST_NAME_OFFSET);
  e->beg_voff = get_u64le(in + BNI_ENTRY_BEG_VOFF_OFFSET);
  e->end_voff = get_u64le(in + BNI_ENTRY_END_VOFF_OFFSET);
  e->n_records = get_u32le(in + BNI_ENTRY_N_RECORDS_OFFSET);
  e->reserved = get_u32le(in + BNI_ENTRY_RESERVED_OFFSET);
}

static int validate_header(const bni_file_header_t *h) {
  if (h->version != BNI_FORMAT_VERSION) {
    bni_print_error("unsupported BNI version %u", h->version);
    return -1;
  }
  if (h->header_size != BNI_HEADER_SIZE) {
    bni_print_error("unsupported BNI header size %u", h->header_size);
    return -1;
  }
  if (h->entry_size != BNI_ENTRY_SIZE) {
    bni_print_error("unsupported BNI entry size %u", h->entry_size);
    return -1;
  }
  if ((h->flags & BNI_FLAG_BGZF_BLOCKS) == 0) {
    bni_print_error("unsupported BNI mode flags 0x%08x", h->flags);
    return -1;
  }
  if (h->entries_offset < BNI_HEADER_SIZE) {
    bni_print_error("invalid entries offset");
    return -1;
  }
  if (h->strings_offset < h->entries_offset) {
    bni_print_error("invalid string-table offset");
    return -1;
  }
  if (h->n_blocks > (UINT64_MAX / BNI_ENTRY_SIZE)) {
    bni_print_error("invalid entry count");
    return -1;
  }
  uint64_t entry_bytes = h->n_blocks * (uint64_t)BNI_ENTRY_SIZE;
  if (h->entries_offset > UINT64_MAX - entry_bytes ||
      h->strings_offset < h->entries_offset + entry_bytes) {
    bni_print_error("invalid entries/string-table layout");
    return -1;
  }
  return 0;
}

static const char *string_at_checked(const bni_index_t *idx, uint64_t off, uint64_t entry_i,
                                     const char *which) {
  if (off >= idx->header.strings_size) {
    bni_print_error("invalid %s name offset in entry %" PRIu64, which, entry_i);
    return NULL;
  }
  const char *name = idx->strings + off;
  size_t max_len = (size_t)(idx->header.strings_size - off);
  if (memchr(name, '\0', max_len) == NULL) {
    bni_print_error("unterminated %s name in entry %" PRIu64, which, entry_i);
    return NULL;
  }
  return name;
}

int bni_write_index_file(const char *path, const bni_file_header_t *header,
                         const bni_entry_t *entries, const char *strings) {
  FILE *fp = fopen(path, "wb");
  if (fp == NULL) {
    bni_print_error("could not open %s for writing: %s", path, strerror(errno));
    return -1;
  }
  unsigned char hbuf[BNI_HEADER_SIZE];
  encode_header(hbuf, header);
  if (write_exact(fp, hbuf, sizeof(hbuf)) != 0) {
    bni_print_error("failed writing BNI header to %s", path);
    close_ignoring_error(fp);
    return -1;
  }
  unsigned char *ebuf = NULL;
  if (header->n_blocks > 0) {
    ebuf = (unsigned char *)malloc((size_t)BNI_ENTRY_WRITE_CHUNK * (size_t)BNI_ENTRY_SIZE);
    if (ebuf == NULL) {
      bni_print_error("out of memory while writing BNI entries");
      close_ignoring_error(fp);
      return -1;
    }
  }
  for (uint64_t i = 0; i < header->n_blocks;) {
    uint64_t remaining = header->n_blocks - i;
    size_t n_entries =
        remaining < BNI_ENTRY_WRITE_CHUNK ? (size_t)remaining : (size_t)BNI_ENTRY_WRITE_CHUNK;
    for (size_t j = 0; j < n_entries; ++j) {
      encode_entry(ebuf + j * (size_t)BNI_ENTRY_SIZE, &entries[i + j]);
    }
    if (write_exact(fp, ebuf, n_entries * (size_t)BNI_ENTRY_SIZE) != 0) {
      bni_print_error("failed writing BNI entries to %s", path);
      free(ebuf);
      close_ignoring_error(fp);
      return -1;
    }
    i += (uint64_t)n_entries;
  }
  free(ebuf);
  if (header->strings_size > 0 && write_exact(fp, strings, (size_t)header->strings_size) != 0) {
    bni_print_error("failed writing BNI string table to %s", path);
    close_ignoring_error(fp);
    return -1;
  }
  if (fclose(fp) != 0) {
    bni_print_error("failed closing %s: %s", path, strerror(errno));
    return -1;
  }
  return 0;
}

int bni_load_index_file(const char *path, bni_index_t *idx) {
  memset(idx, 0, sizeof(*idx));
  FILE *fp = fopen(path, "rb");
  if (fp == NULL) {
    bni_print_error("could not open %s: %s", path, strerror(errno));
    return -1;
  }
  unsigned char hbuf[BNI_HEADER_SIZE];
  if (read_exact(fp, hbuf, sizeof(hbuf)) != 0) {
    bni_print_error("failed reading BNI header from %s", path);
    close_ignoring_error(fp);
    return -1;
  }
  if (decode_header(hbuf, &idx->header) != 0 || validate_header(&idx->header) != 0) {
    bni_print_error("%s is not a valid BNI v2 BGZF-block index", path);
    close_ignoring_error(fp);
    return -1;
  }
  uint64_t entry_bytes_u64 = idx->header.n_blocks * (uint64_t)BNI_ENTRY_SIZE;
  if (entry_bytes_u64 > SIZE_MAX || idx->header.strings_size > SIZE_MAX - 1) {
    bni_print_error("index is too large for this platform");
    close_ignoring_error(fp);
    return -1;
  }
  if (idx->header.n_blocks > 0) {
    idx->entries = (bni_entry_t *)calloc((size_t)idx->header.n_blocks, sizeof(bni_entry_t));
    if (idx->entries == NULL) {
      bni_print_error("out of memory while loading entries");
      close_ignoring_error(fp);
      return -1;
    }
  }
  if (fseeko(fp, (off_t)idx->header.entries_offset, SEEK_SET) != 0) {
    bni_print_error("failed seeking to entries in %s", path);
    bni_index_destroy(idx);
    close_ignoring_error(fp);
    return -1;
  }
  unsigned char ebuf[BNI_ENTRY_SIZE];
  for (uint64_t i = 0; i < idx->header.n_blocks; ++i) {
    if (read_exact(fp, ebuf, sizeof(ebuf)) != 0) {
      bni_print_error("failed reading entry from %s", path);
      bni_index_destroy(idx);
      close_ignoring_error(fp);
      return -1;
    }
    decode_entry(ebuf, &idx->entries[i]);
  }
  idx->strings = (char *)calloc((size_t)idx->header.strings_size + 1, 1);
  if (idx->strings == NULL) {
    bni_print_error("out of memory while loading string table");
    bni_index_destroy(idx);
    close_ignoring_error(fp);
    return -1;
  }
  if (idx->header.strings_size > 0) {
    if (fseeko(fp, (off_t)idx->header.strings_offset, SEEK_SET) != 0) {
      bni_print_error("failed seeking to string table in %s", path);
      bni_index_destroy(idx);
      close_ignoring_error(fp);
      return -1;
    }
    if (read_exact(fp, idx->strings, (size_t)idx->header.strings_size) != 0) {
      bni_print_error("failed reading string table from %s", path);
      bni_index_destroy(idx);
      close_ignoring_error(fp);
      return -1;
    }
  }
  for (uint64_t i = 0; i < idx->header.n_blocks; ++i) {
    const bni_entry_t *e = &idx->entries[i];
    const char *first = string_at_checked(idx, e->first_name_offset, i, "first");
    const char *last = string_at_checked(idx, e->last_name_offset, i, "last");
    if (first == NULL || last == NULL) {
      bni_index_destroy(idx);
      close_ignoring_error(fp);
      return -1;
    }
    if (strcmp(first, last) > 0) {
      bni_print_error("entry %" PRIu64 " has first QNAME greater than last QNAME", i);
      bni_index_destroy(idx);
      close_ignoring_error(fp);
      return -1;
    }
    if (e->beg_voff >= e->end_voff) {
      bni_print_error("entry %" PRIu64 " has an empty virtual-offset range", i);
      bni_index_destroy(idx);
      close_ignoring_error(fp);
      return -1;
    }
    if (e->n_records == 0) {
      bni_print_error("entry %" PRIu64 " has zero records", i);
      bni_index_destroy(idx);
      close_ignoring_error(fp);
      return -1;
    }
    if (i > 0) {
      const bni_entry_t *prev_e = &idx->entries[i - 1];
      const char *prev_first = bni_entry_first_name(idx, prev_e);
      const char *prev_last = bni_entry_last_name(idx, prev_e);
      if (strcmp(prev_first, first) > 0) {
        bni_print_error("entries are not sorted by first QNAME");
        bni_index_destroy(idx);
        close_ignoring_error(fp);
        return -1;
      }
      if (strcmp(prev_last, last) > 0) {
        bni_print_error("entries are not sorted by last QNAME");
        bni_index_destroy(idx);
        close_ignoring_error(fp);
        return -1;
      }
      if (prev_e->end_voff > e->beg_voff) {
        bni_print_error("entries have overlapping virtual-offset ranges");
        bni_index_destroy(idx);
        close_ignoring_error(fp);
        return -1;
      }
    }
  }
  if (fclose(fp) != 0) {
    bni_print_error("failed closing %s: %s", path, strerror(errno));
    bni_index_destroy(idx);
    return -1;
  }
  return 0;
}

void bni_index_destroy(bni_index_t *idx) {
  if (idx == NULL) {
    return;
  }
  free(idx->entries);
  free(idx->strings);
  memset(idx, 0, sizeof(*idx));
}

const char *bni_entry_first_name(const bni_index_t *idx, const bni_entry_t *entry) {
  if (idx == NULL || entry == NULL || idx->strings == NULL) {
    return NULL;
  }
  if (entry->first_name_offset >= idx->header.strings_size) {
    return NULL;
  }
  return idx->strings + entry->first_name_offset;
}

const char *bni_entry_last_name(const bni_index_t *idx, const bni_entry_t *entry) {
  if (idx == NULL || entry == NULL || idx->strings == NULL) {
    return NULL;
  }
  if (entry->last_name_offset >= idx->header.strings_size) {
    return NULL;
  }
  return idx->strings + entry->last_name_offset;
}

bni_index_t *bni_index_open(const char *path) {
  bni_index_t *idx = (bni_index_t *)calloc(1, sizeof(*idx));
  if (idx == NULL) {
    bni_print_error("out of memory while allocating index");
    return NULL;
  }
  if (bni_load_index_file(path, idx) != 0) {
    free(idx);
    return NULL;
  }
  return idx;
}

void bni_index_close(bni_index_t *idx) {
  if (idx == NULL) {
    return;
  }
  bni_index_destroy(idx);
  free(idx);
}

const bni_entry_t *bni_find_entry(const bni_index_t *idx, const char *name) {
  if (idx == NULL || name == NULL) {
    return NULL;
  }
  uint64_t lo = 0;
  uint64_t hi = idx->header.n_blocks;
  while (lo < hi) {
    uint64_t mid = lo + (hi - lo) / 2;
    const char *mid_last = bni_entry_last_name(idx, &idx->entries[mid]);
    int cmp = strcmp(mid_last, name);
    if (cmp < 0) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  if (lo == idx->header.n_blocks) {
    return NULL;
  }
  return &idx->entries[lo];
}
