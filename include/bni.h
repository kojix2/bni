#ifndef BNI_H
#define BNI_H

#include <stdint.h>
#include <stddef.h>

#define BNI_VERSION_STRING "0.0.0"

#define BNI_MAGIC0 'B'
#define BNI_MAGIC1 'N'
#define BNI_MAGIC2 'I'
#define BNI_MAGIC3 1

#define BNI_FORMAT_VERSION 1u
#define BNI_HEADER_SIZE 128u
#define BNI_ENTRY_SIZE 32u

#define BNI_SORT_QUERYNAME_LEX 1u
#define BNI_FLAG_QUERYNAME_GROUPED 0x00000001u

typedef struct {
    uint32_t version;
    uint32_t header_size;
    uint32_t flags;
    uint64_t n_names;
    uint64_t n_records;
    uint64_t entries_offset;
    uint64_t strings_offset;
    uint64_t strings_size;
    uint64_t bam_size;
    int64_t  bam_mtime;
    uint64_t header_hash;
    uint32_t sort_order;
    uint32_t entry_size;
} bni_file_header_t;

typedef struct {
    uint64_t name_offset;
    uint64_t beg_voff;
    uint64_t end_voff;
    uint32_t n_records;
    uint32_t reserved;
} bni_entry_t;

typedef struct {
    bni_file_header_t header;
    bni_entry_t *entries;
    char *strings;
} bni_index_t;

int bni_write_index_file(const char *path, const bni_file_header_t *header,
                         const bni_entry_t *entries, const char *strings);
int bni_load_index_file(const char *path, bni_index_t *idx);
void bni_index_destroy(bni_index_t *idx);
const char *bni_entry_name(const bni_index_t *idx, const bni_entry_t *entry);
const bni_entry_t *bni_find_entry(const bni_index_t *idx, const char *name);

char *bni_default_index_path(const char *bam_path);
uint64_t bni_fnv1a64(const void *data, size_t len);
int bni_file_metadata(const char *path, uint64_t *size_out, int64_t *mtime_out);
int bni_path_exists(const char *path);
int bni_has_suffix(const char *s, const char *suffix);
int bni_parse_threads(const char *s, int *out);
char *bni_strdup(const char *s);
void bni_print_error(const char *fmt, ...);
void bni_print_warning(const char *fmt, ...);
void bni_format_u64(char *buf, size_t buflen, uint64_t value);

int bni_cmd_index(int argc, char **argv);
int bni_cmd_get(int argc, char **argv);
int bni_cmd_stats(int argc, char **argv);
int bni_cmd_check(int argc, char **argv);

#endif
