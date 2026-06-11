#define _POSIX_C_SOURCE 200809L

#include "bni_internal.h"

#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

static void put_u32le(unsigned char *p, uint32_t v) {
    p[0] = (unsigned char)(v & 0xffu);
    p[1] = (unsigned char)((v >> 8) & 0xffu);
    p[2] = (unsigned char)((v >> 16) & 0xffu);
    p[3] = (unsigned char)((v >> 24) & 0xffu);
}

static void put_u64le(unsigned char *p, uint64_t v) {
    for (int i = 0; i < 8; ++i) p[i] = (unsigned char)((v >> (8 * i)) & 0xffu);
}

static uint32_t get_u32le(const unsigned char *p) {
    return ((uint32_t)p[0]) | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}

static uint64_t get_u64le(const unsigned char *p) {
    uint64_t v = 0;
    for (int i = 7; i >= 0; --i) { v <<= 8; v |= (uint64_t)p[i]; }
    return v;
}

static int write_exact(FILE *fp, const void *buf, size_t n) {
    return fwrite(buf, 1, n, fp) == n ? 0 : -1;
}

static int read_exact(FILE *fp, void *buf, size_t n) {
    return fread(buf, 1, n, fp) == n ? 0 : -1;
}

static void encode_header(unsigned char out[BNI_HEADER_SIZE], const bni_file_header_t *h) {
    memset(out, 0, BNI_HEADER_SIZE);
    out[0] = BNI_MAGIC0;
    out[1] = BNI_MAGIC1;
    out[2] = BNI_MAGIC2;
    out[3] = BNI_MAGIC3;
    put_u32le(out + 4, h->version);
    put_u32le(out + 8, h->header_size);
    put_u32le(out + 12, h->flags);
    put_u64le(out + 16, h->n_names);
    put_u64le(out + 24, h->n_records);
    put_u64le(out + 32, h->entries_offset);
    put_u64le(out + 40, h->strings_offset);
    put_u64le(out + 48, h->strings_size);
    put_u64le(out + 56, h->bam_size);
    put_u64le(out + 64, (uint64_t)h->bam_mtime);
    put_u64le(out + 72, h->header_hash);
    put_u32le(out + 80, h->sort_order);
    put_u32le(out + 84, h->entry_size);
}

static int decode_header(const unsigned char in[BNI_HEADER_SIZE], bni_file_header_t *h) {
    if (in[0] != BNI_MAGIC0 || in[1] != BNI_MAGIC1 || in[2] != BNI_MAGIC2 || in[3] != BNI_MAGIC3) return -1;
    memset(h, 0, sizeof(*h));
    h->version = get_u32le(in + 4);
    h->header_size = get_u32le(in + 8);
    h->flags = get_u32le(in + 12);
    h->n_names = get_u64le(in + 16);
    h->n_records = get_u64le(in + 24);
    h->entries_offset = get_u64le(in + 32);
    h->strings_offset = get_u64le(in + 40);
    h->strings_size = get_u64le(in + 48);
    h->bam_size = get_u64le(in + 56);
    h->bam_mtime = (int64_t)get_u64le(in + 64);
    h->header_hash = get_u64le(in + 72);
    h->sort_order = get_u32le(in + 80);
    h->entry_size = get_u32le(in + 84);
    return 0;
}

static void encode_entry(unsigned char out[BNI_ENTRY_SIZE], const bni_entry_t *e) {
    memset(out, 0, BNI_ENTRY_SIZE);
    put_u64le(out + 0, e->name_offset);
    put_u64le(out + 8, e->beg_voff);
    put_u64le(out + 16, e->end_voff);
    put_u32le(out + 24, e->n_records);
    put_u32le(out + 28, e->reserved);
}

static void decode_entry(const unsigned char in[BNI_ENTRY_SIZE], bni_entry_t *e) {
    e->name_offset = get_u64le(in + 0);
    e->beg_voff = get_u64le(in + 8);
    e->end_voff = get_u64le(in + 16);
    e->n_records = get_u32le(in + 24);
    e->reserved = get_u32le(in + 28);
}

static int validate_header(const bni_file_header_t *h) {
    if (h->version != BNI_FORMAT_VERSION) { bni_print_error("unsupported BNI version %u", h->version); return -1; }
    if (h->header_size != BNI_HEADER_SIZE) { bni_print_error("unsupported BNI header size %u", h->header_size); return -1; }
    if (h->entry_size != BNI_ENTRY_SIZE) { bni_print_error("unsupported BNI entry size %u", h->entry_size); return -1; }
    if (h->entries_offset < BNI_HEADER_SIZE) { bni_print_error("invalid entries offset"); return -1; }
    if (h->strings_offset < h->entries_offset) { bni_print_error("invalid string-table offset"); return -1; }
    if (h->n_names > (UINT64_MAX / BNI_ENTRY_SIZE)) { bni_print_error("invalid entry count"); return -1; }
    return 0;
}

int bni_write_index_file(const char *path, const bni_file_header_t *header,
                         const bni_entry_t *entries, const char *strings) {
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) { bni_print_error("could not open %s for writing: %s", path, strerror(errno)); return -1; }
    unsigned char hbuf[BNI_HEADER_SIZE];
    encode_header(hbuf, header);
    if (write_exact(fp, hbuf, sizeof(hbuf)) != 0) { bni_print_error("failed writing BNI header to %s", path); fclose(fp); return -1; }
    unsigned char ebuf[BNI_ENTRY_SIZE];
    for (uint64_t i = 0; i < header->n_names; ++i) {
        encode_entry(ebuf, &entries[i]);
        if (write_exact(fp, ebuf, sizeof(ebuf)) != 0) { bni_print_error("failed writing BNI entries to %s", path); fclose(fp); return -1; }
    }
    if (header->strings_size > 0 && write_exact(fp, strings, (size_t)header->strings_size) != 0) {
        bni_print_error("failed writing BNI string table to %s", path); fclose(fp); return -1;
    }
    if (fclose(fp) != 0) { bni_print_error("failed closing %s: %s", path, strerror(errno)); return -1; }
    return 0;
}

int bni_load_index_file(const char *path, bni_index_t *idx) {
    memset(idx, 0, sizeof(*idx));
    FILE *fp = fopen(path, "rb");
    if (fp == NULL) { bni_print_error("could not open %s: %s", path, strerror(errno)); return -1; }
    unsigned char hbuf[BNI_HEADER_SIZE];
    if (read_exact(fp, hbuf, sizeof(hbuf)) != 0) { bni_print_error("failed reading BNI header from %s", path); fclose(fp); return -1; }
    if (decode_header(hbuf, &idx->header) != 0 || validate_header(&idx->header) != 0) {
        bni_print_error("%s is not a valid BNI v1 index", path); fclose(fp); return -1;
    }
    uint64_t entry_bytes_u64 = idx->header.n_names * (uint64_t)BNI_ENTRY_SIZE;
    if (entry_bytes_u64 > SIZE_MAX || idx->header.strings_size > SIZE_MAX - 1) {
        bni_print_error("index is too large for this platform"); fclose(fp); return -1;
    }
    if (idx->header.n_names > 0) {
        idx->entries = (bni_entry_t *)calloc((size_t)idx->header.n_names, sizeof(bni_entry_t));
        if (idx->entries == NULL) { bni_print_error("out of memory while loading entries"); fclose(fp); return -1; }
    }
    if (fseeko(fp, (off_t)idx->header.entries_offset, SEEK_SET) != 0) {
        bni_print_error("failed seeking to entries in %s", path); bni_index_destroy(idx); fclose(fp); return -1;
    }
    unsigned char ebuf[BNI_ENTRY_SIZE];
    for (uint64_t i = 0; i < idx->header.n_names; ++i) {
        if (read_exact(fp, ebuf, sizeof(ebuf)) != 0) { bni_print_error("failed reading entry from %s", path); bni_index_destroy(idx); fclose(fp); return -1; }
        decode_entry(ebuf, &idx->entries[i]);
    }
    idx->strings = (char *)calloc((size_t)idx->header.strings_size + 1, 1);
    if (idx->strings == NULL) { bni_print_error("out of memory while loading string table"); bni_index_destroy(idx); fclose(fp); return -1; }
    if (idx->header.strings_size > 0) {
        if (fseeko(fp, (off_t)idx->header.strings_offset, SEEK_SET) != 0) {
            bni_print_error("failed seeking to string table in %s", path); bni_index_destroy(idx); fclose(fp); return -1;
        }
        if (read_exact(fp, idx->strings, (size_t)idx->header.strings_size) != 0) {
            bni_print_error("failed reading string table from %s", path); bni_index_destroy(idx); fclose(fp); return -1;
        }
    }
    for (uint64_t i = 0; i < idx->header.n_names; ++i) {
        const bni_entry_t *e = &idx->entries[i];
        if (e->name_offset >= idx->header.strings_size) { bni_print_error("invalid name offset in entry %" PRIu64, i); bni_index_destroy(idx); fclose(fp); return -1; }
        const char *name = idx->strings + e->name_offset;
        size_t max_len = (size_t)(idx->header.strings_size - e->name_offset);
        if (memchr(name, '\0', max_len) == NULL) { bni_print_error("unterminated name in entry %" PRIu64, i); bni_index_destroy(idx); fclose(fp); return -1; }
        if (i > 0) {
            const char *prev = bni_entry_name(idx, &idx->entries[i - 1]);
            if (strcmp(prev, name) >= 0) { bni_print_error("entries are not strictly sorted by QNAME"); bni_index_destroy(idx); fclose(fp); return -1; }
        }
    }
    fclose(fp);
    return 0;
}

void bni_index_destroy(bni_index_t *idx) {
    if (idx == NULL) return;
    free(idx->entries);
    free(idx->strings);
    memset(idx, 0, sizeof(*idx));
}

const char *bni_entry_name(const bni_index_t *idx, const bni_entry_t *entry) {
    if (idx == NULL || entry == NULL || idx->strings == NULL) return NULL;
    if (entry->name_offset >= idx->header.strings_size) return NULL;
    return idx->strings + entry->name_offset;
}

const bni_entry_t *bni_find_entry(const bni_index_t *idx, const char *name) {
    if (idx == NULL || name == NULL) return NULL;
    uint64_t lo = 0, hi = idx->header.n_names;
    while (lo < hi) {
        uint64_t mid = lo + (hi - lo) / 2;
        const char *mid_name = bni_entry_name(idx, &idx->entries[mid]);
        int cmp = strcmp(mid_name, name);
        if (cmp < 0) lo = mid + 1;
        else if (cmp > 0) hi = mid;
        else return &idx->entries[mid];
    }
    return NULL;
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
    if (idx == NULL) return;
    bni_index_destroy(idx);
    free(idx);
}
