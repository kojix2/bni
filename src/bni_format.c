#include "bni_internal.h"

#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef _WIN32
#include <sys/mman.h>
#include <unistd.h>
#endif
#include <sys/types.h>

#define BNI_ENTRY_WRITE_CHUNK 8192u

enum {
  BNI_DEFAULT_MMAP_PAGE_SIZE = 4096,
  BNI_TEMP_CREATE_ATTEMPTS = 100,
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
  BNI_HEADER_BAM_MTIME_NSEC_OFFSET = 88,

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
  put_u32le(out + BNI_HEADER_BAM_MTIME_NSEC_OFFSET, h->bam_mtime_nsec);
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
  h->bam_mtime_nsec = get_u32le(in + BNI_HEADER_BAM_MTIME_NSEC_OFFSET);
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
  const uint32_t known_flags = BNI_FLAG_BGZF_BLOCKS | BNI_FLAG_MTIME_NSEC;
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
  if ((h->flags & ~known_flags) != 0) {
    bni_print_error("unknown BNI flag bits 0x%08x", h->flags & ~known_flags);
    return -1;
  }
  if (h->sort_order != BNI_SORT_QUERYNAME_LEX) {
    bni_print_error("unsupported BNI sort order %u", h->sort_order);
    return -1;
  }
  if ((h->flags & BNI_FLAG_MTIME_NSEC) != 0 && h->bam_mtime_nsec >= 1000000000u) {
    bni_print_error("invalid BAM mtime nanoseconds %u", h->bam_mtime_nsec);
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

static int validate_string_table_bounds(const bni_index_t *idx) {
  if (idx->header.n_blocks == 0) {
    return 0;
  }
  if (idx->header.strings_size == 0 || idx->strings == NULL) {
    bni_print_error("index entries require a non-empty string table");
    return -1;
  }
  // NOLINTNEXTLINE(clang-analyzer-security.ArrayBound)
  if (idx->strings[idx->header.strings_size - 1] != '\0') {
    bni_print_error("index string table is not NUL-terminated");
    return -1;
  }
  return 0;
}

static int write_index_contents(FILE *fp, const char *path, const bni_file_header_t *header,
                                const bni_entry_t *entries, const char *strings) {
  unsigned char hbuf[BNI_HEADER_SIZE];
  encode_header(hbuf, header);
  if (write_exact(fp, hbuf, sizeof(hbuf)) != 0) {
    bni_print_error("failed writing BNI header to %s", path);
    return -1;
  }
  unsigned char *ebuf = NULL;
  if (header->n_blocks > 0) {
    ebuf = (unsigned char *)malloc((size_t)BNI_ENTRY_WRITE_CHUNK * (size_t)BNI_ENTRY_SIZE);
    if (ebuf == NULL) {
      bni_print_error("out of memory while writing BNI entries");
      return -1;
    }
  }
  for (uint64_t i = 0; i < header->n_blocks;) {
    uint64_t remaining = header->n_blocks - i;
    size_t n_entries =
        remaining < BNI_ENTRY_WRITE_CHUNK ? (size_t)remaining : (size_t)BNI_ENTRY_WRITE_CHUNK;
    for (size_t j = 0; j < n_entries; ++j) {
      encode_entry(ebuf + (j * (size_t)BNI_ENTRY_SIZE), &entries[i + j]);
    }
    if (write_exact(fp, ebuf, n_entries * (size_t)BNI_ENTRY_SIZE) != 0) {
      bni_print_error("failed writing BNI entries to %s", path);
      free(ebuf);
      return -1;
    }
    i += (uint64_t)n_entries;
  }
  free(ebuf);
  if (header->strings_size > 0 && write_exact(fp, strings, (size_t)header->strings_size) != 0) {
    bni_print_error("failed writing BNI string table to %s", path);
    return -1;
  }
  return 0;
}

static int rewind_spool(FILE *spool, const char *description) {
  if (fseeko(spool, 0, SEEK_SET) != 0) {
    bni_print_error("failed rewinding temporary index %s: %s", description, strerror(errno));
    return -1;
  }
  return 0;
}

static int write_spooled_index_contents(FILE *fp, const char *path,
                                        const bni_file_header_t *header, FILE *entries_spool,
                                        FILE *strings_spool) {
  unsigned char hbuf[BNI_HEADER_SIZE];
  encode_header(hbuf, header);
  if (write_exact(fp, hbuf, sizeof(hbuf)) != 0) {
    bni_print_error("failed writing BNI header to %s", path);
    return -1;
  }
  if (rewind_spool(entries_spool, "entry table") != 0) {
    return -1;
  }
  unsigned char ebuf[BNI_ENTRY_SIZE];
  for (uint64_t i = 0; i < header->n_blocks; ++i) {
    bni_entry_t entry;
    if (fread(&entry, sizeof(entry), 1, entries_spool) != 1) {
      bni_print_error("failed reading temporary index entry table");
      return -1;
    }
    encode_entry(ebuf, &entry);
    if (write_exact(fp, ebuf, sizeof(ebuf)) != 0) {
      bni_print_error("failed writing BNI entries to %s", path);
      return -1;
    }
  }
  if (rewind_spool(strings_spool, "string table") != 0) {
    return -1;
  }
  unsigned char buffer[64 * 1024];
  uint64_t remaining = header->strings_size;
  while (remaining > 0) {
    size_t chunk = remaining < sizeof(buffer) ? (size_t)remaining : sizeof(buffer);
    if (fread(buffer, 1, chunk, strings_spool) != chunk) {
      bni_print_error("failed reading temporary index string table");
      return -1;
    }
    if (write_exact(fp, buffer, chunk) != 0) {
      bni_print_error("failed writing BNI string table to %s", path);
      return -1;
    }
    remaining -= (uint64_t)chunk;
  }
  return 0;
}

static int create_index_temp(const char *path, char **temp_path_out) {
  static const char suffix_format[] = ".tmp.%ld.%u";
  size_t path_len = strlen(path);
  if (path_len > SIZE_MAX - sizeof(suffix_format) - 32) {
    errno = ENAMETOOLONG;
    return -1;
  }
  size_t temp_path_size = path_len + sizeof(suffix_format) + 32;
  char *temp_path = (char *)malloc(temp_path_size);
  if (temp_path == NULL) {
    errno = ENOMEM;
    return -1;
  }

  for (unsigned int attempt = 0; attempt < BNI_TEMP_CREATE_ATTEMPTS; ++attempt) {
    (void)snprintf(temp_path, temp_path_size, "%s.tmp.%ld.%u", path, (long)getpid(), attempt);
    int fd = open(temp_path, O_WRONLY | O_CREAT | O_EXCL, 0666);
    if (fd >= 0) {
      *temp_path_out = temp_path;
      return fd;
    }
    if (errno != EEXIST) {
      break;
    }
  }
  int saved_errno = errno;
  free(temp_path);
  errno = saved_errno;
  return -1;
}

static int publish_index_temp(const char *temp_path, const char *path, int replace) {
  if (replace) {
    return rename(temp_path, path);
  }
  if (link(temp_path, path) != 0) {
    return -1;
  }
  if (unlink(temp_path) != 0) {
    return -1;
  }
  return 0;
}

static int write_index_file_atomic(const char *path, const bni_file_header_t *header,
                                   const bni_entry_t *entries, const char *strings, int replace) {
  char *temp_path = NULL;
  int fd = create_index_temp(path, &temp_path);
  if (fd < 0) {
    bni_print_error("could not create temporary index for %s: %s", path, strerror(errno));
    return -1;
  }
  FILE *fp = fdopen(fd, "wb");
  if (fp == NULL) {
    int saved_errno = errno;
    (void)close(fd);
    (void)unlink(temp_path);
    bni_print_error("could not open temporary index for %s: %s", path, strerror(saved_errno));
    free(temp_path);
    return -1;
  }

  int status = write_index_contents(fp, path, header, entries, strings);
  if (status == 0 && fflush(fp) != 0) {
    bni_print_error("failed flushing temporary index for %s: %s", path, strerror(errno));
    status = -1;
  }
  if (status == 0 && fsync(fileno(fp)) != 0) {
    bni_print_error("failed syncing temporary index for %s: %s", path, strerror(errno));
    status = -1;
  }
  if (fclose(fp) != 0) {
    if (status == 0) {
      bni_print_error("failed closing temporary index for %s: %s", path, strerror(errno));
    }
    status = -1;
  }
  if (status == 0 && publish_index_temp(temp_path, path, replace) != 0) {
    bni_print_error("failed publishing index %s: %s", path, strerror(errno));
    status = -1;
  }
  if (status != 0) {
    (void)unlink(temp_path);
  }
  free(temp_path);
  return status;
}

int bni_write_spooled_index_file(const char *path, const bni_file_header_t *header,
                                 FILE *entries_spool, FILE *strings_spool, int replace) {
  char *temp_path = NULL;
  int fd = create_index_temp(path, &temp_path);
  if (fd < 0) {
    bni_print_error("could not create temporary index for %s: %s", path, strerror(errno));
    return -1;
  }
  FILE *fp = fdopen(fd, "wb");
  if (fp == NULL) {
    int saved_errno = errno;
    (void)close(fd);
    (void)unlink(temp_path);
    bni_print_error("could not open temporary index for %s: %s", path, strerror(saved_errno));
    free(temp_path);
    return -1;
  }

  int status = write_spooled_index_contents(fp, path, header, entries_spool, strings_spool);
  if (status == 0 && fflush(fp) != 0) {
    bni_print_error("failed flushing temporary index for %s: %s", path, strerror(errno));
    status = -1;
  }
  if (status == 0 && fsync(fileno(fp)) != 0) {
    bni_print_error("failed syncing temporary index for %s: %s", path, strerror(errno));
    status = -1;
  }
  if (fclose(fp) != 0) {
    if (status == 0) {
      bni_print_error("failed closing temporary index for %s: %s", path, strerror(errno));
    }
    status = -1;
  }
  if (status == 0 && publish_index_temp(temp_path, path, replace) != 0) {
    bni_print_error("failed publishing index %s: %s", path, strerror(errno));
    status = -1;
  }
  if (status != 0) {
    (void)unlink(temp_path);
  }
  free(temp_path);
  return status;
}

int bni_write_index_file(const char *path, const bni_file_header_t *header,
                         const bni_entry_t *entries, const char *strings) {
  return write_index_file_atomic(path, header, entries, strings, 1);
}

int bni_write_index_file_exclusive(const char *path, const bni_file_header_t *header,
                                   const bni_entry_t *entries, const char *strings) {
  return write_index_file_atomic(path, header, entries, strings, 0);
}

static int read_index_header(FILE *fp, const char *path, bni_file_header_t *header) {
  unsigned char hbuf[BNI_HEADER_SIZE];
  if (read_exact(fp, hbuf, sizeof(hbuf)) != 0) {
    bni_print_error("failed reading BNI header from %s", path);
    return -1;
  }
  if (decode_header(hbuf, header) != 0 || validate_header(header) != 0) {
    bni_print_error("%s is not a valid BNI v1 BGZF-block index", path);
    return -1;
  }
  uint64_t entry_bytes_u64 = header->n_blocks * (uint64_t)BNI_ENTRY_SIZE;
  if (entry_bytes_u64 > SIZE_MAX || header->strings_size > SIZE_MAX - 1) {
    bni_print_error("index is too large for this platform");
    return -1;
  }
  return 0;
}

static int validate_index_file_size(FILE *fp, const char *path, const bni_file_header_t *header) {
  if (fseeko(fp, 0, SEEK_END) != 0) {
    bni_print_error("failed seeking to end of %s", path);
    return -1;
  }
  off_t end_pos = ftello(fp);
  if (end_pos < 0) {
    bni_print_error("failed determining size of %s", path);
    return -1;
  }
  uint64_t file_size = (uint64_t)end_pos;
  if (header->strings_offset > UINT64_MAX - header->strings_size ||
      header->strings_offset + header->strings_size > file_size) {
    bni_print_error("index string table extends past end of %s", path);
    return -1;
  }
  return 0;
}

static int read_index_entries(FILE *fp, const char *path, bni_index_t *idx) {
  if (idx->header.n_blocks > 0) {
    idx->entries = (bni_entry_t *)calloc((size_t)idx->header.n_blocks, sizeof(bni_entry_t));
    if (idx->entries == NULL) {
      bni_print_error("out of memory while loading entries");
      return -1;
    }
  }
  if (fseeko(fp, (off_t)idx->header.entries_offset, SEEK_SET) != 0) {
    bni_print_error("failed seeking to entries in %s", path);
    return -1;
  }
  unsigned char ebuf[BNI_ENTRY_SIZE];
  for (uint64_t i = 0; i < idx->header.n_blocks; ++i) {
    if (read_exact(fp, ebuf, sizeof(ebuf)) != 0) {
      bni_print_error("failed reading entry from %s", path);
      return -1;
    }
    decode_entry(ebuf, &idx->entries[i]);
  }
  return 0;
}

static int read_index_strings(FILE *fp, const char *path, bni_index_t *idx) {
  if (idx->header.strings_size > SIZE_MAX - 1) {
    bni_print_error("index string table is too large for this platform");
    return -1;
  }
  idx->strings = (char *)calloc((size_t)idx->header.strings_size + 1, 1);
  if (idx->strings == NULL) {
    bni_print_error("out of memory while loading string table");
    return -1;
  }
  idx->owns_strings = 1;
  if (idx->header.strings_size == 0) {
    return 0;
  }
  if (fseeko(fp, (off_t)idx->header.strings_offset, SEEK_SET) != 0) {
    bni_print_error("failed seeking to string table in %s", path);
    return -1;
  }
  if (read_exact(fp, idx->strings, (size_t)idx->header.strings_size) != 0) {
    bni_print_error("failed reading string table from %s", path);
    return -1;
  }
  return 0;
}

static int map_index_strings(FILE *fp, bni_index_t *idx) {
#ifdef _WIN32
  (void)fp;
  (void)idx;
  return -1;
#else
  if (idx->header.strings_size == 0) {
    return 0;
  }
  long page_size_long = sysconf(_SC_PAGE_SIZE);
  size_t page_size =
      page_size_long > 0 ? (size_t)page_size_long : (size_t)BNI_DEFAULT_MMAP_PAGE_SIZE;
  uint64_t map_offset = idx->header.strings_offset - (idx->header.strings_offset % page_size);
  uint64_t map_delta = idx->header.strings_offset - map_offset;
  if (map_delta > SIZE_MAX || idx->header.strings_size > SIZE_MAX - (size_t)map_delta) {
    return -1;
  }
  size_t map_size = (size_t)map_delta + (size_t)idx->header.strings_size;
  if (map_offset > (uint64_t)INT64_MAX) {
    return -1;
  }
  void *mapping = mmap(NULL, map_size, PROT_READ, MAP_PRIVATE, fileno(fp), (off_t)map_offset);
  if (mapping == MAP_FAILED) {
    return -1;
  }
  idx->mapping = mapping;
  idx->mapping_size = map_size;
  idx->strings = (char *)mapping + (size_t)map_delta;
  idx->owns_strings = 0;
  return 0;
#endif
}

static int load_index_strings(FILE *fp, const char *path, bni_index_t *idx) {
  if (map_index_strings(fp, idx) == 0) {
    return 0;
  }
  return read_index_strings(fp, path, idx);
}

static int validate_name_offset(const bni_index_t *idx, uint64_t offset, const char *kind,
                                uint64_t entry_i, const char **name_out) {
  if (offset >= idx->header.strings_size) {
    bni_print_error("invalid %s name offset in entry %" PRIu64, kind, entry_i);
    return -1;
  }
  if (offset > 0 && idx->strings[offset - 1] != '\0') {
    bni_print_error("%s name offset in entry %" PRIu64 " is not at a string boundary", kind,
                    entry_i);
    return -1;
  }
  const char *name = idx->strings + offset;
  size_t remaining = (size_t)(idx->header.strings_size - offset);
  if (memchr(name, '\0', remaining) == NULL) {
    bni_print_error("%s name in entry %" PRIu64 " is not NUL-terminated", kind, entry_i);
    return -1;
  }
  *name_out = name;
  return 0;
}

static int validate_index_entry(const bni_index_t *idx, uint64_t entry_i, const char **first_out,
                                const char **last_out) {
  const bni_entry_t *e = &idx->entries[entry_i];
  if (validate_name_offset(idx, e->first_name_offset, "first", entry_i, first_out) != 0 ||
      validate_name_offset(idx, e->last_name_offset, "last", entry_i, last_out) != 0) {
    return -1;
  }
  if (strcmp(*first_out, *last_out) > 0) {
    bni_print_error("entry %" PRIu64 " has decreasing QNAME range '%s'..'%s'", entry_i,
                    *first_out, *last_out);
    return -1;
  }
  if (e->beg_voff >= e->end_voff) {
    bni_print_error("entry %" PRIu64 " has an empty virtual-offset range", entry_i);
    return -1;
  }
  if (e->n_records == 0) {
    bni_print_error("entry %" PRIu64 " has zero records", entry_i);
    return -1;
  }
  if (entry_i == 0) {
    return 0;
  }

  const bni_entry_t *prev_e = &idx->entries[entry_i - 1];
  if (prev_e->end_voff != e->beg_voff) {
    bni_print_error("entries have non-contiguous virtual-offset ranges");
    return -1;
  }
  return 0;
}

static int validate_index_entries(const bni_index_t *idx) {
  if (validate_string_table_bounds(idx) != 0) {
    return -1;
  }
  uint64_t total_records = 0;
  const char *previous_last = NULL;
  for (uint64_t i = 0; i < idx->header.n_blocks; ++i) {
    const char *first = NULL;
    const char *last = NULL;
    if (validate_index_entry(idx, i, &first, &last) != 0) {
      return -1;
    }
    if (previous_last != NULL && strcmp(previous_last, first) > 0) {
      bni_print_error("QNAME ranges decrease between entries %" PRIu64 " and %" PRIu64
                      " near '%s' -> '%s'",
                      i - 1, i, previous_last, first);
      return -1;
    }
    previous_last = last;
    if (total_records > UINT64_MAX - idx->entries[i].n_records) {
      bni_print_error("index record count overflows");
      return -1;
    }
    total_records += idx->entries[i].n_records;
  }
  if (total_records != idx->header.n_records) {
    bni_print_error("index record count mismatch: header=%" PRIu64 " entries=%" PRIu64,
                    idx->header.n_records, total_records);
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
  if (read_index_header(fp, path, &idx->header) != 0 ||
      validate_index_file_size(fp, path, &idx->header) != 0 ||
      read_index_entries(fp, path, idx) != 0 || load_index_strings(fp, path, idx) != 0 ||
      validate_index_entries(idx) != 0) {
    bni_index_destroy(idx);
    close_ignoring_error(fp);
    return -1;
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
  if (idx->mapping != NULL) {
#ifndef _WIN32
    munmap(idx->mapping, idx->mapping_size);
#endif
  } else if (idx->owns_strings) {
    free(idx->strings);
  }
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
    uint64_t mid = lo + ((hi - lo) / 2);
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
