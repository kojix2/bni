#ifndef BNI_H
#define BNI_H

#include <stdint.h>
#include <stddef.h>

#include <htslib/sam.h>

#if defined(_WIN32) && defined(BNI_BUILD_SHARED)
#  define BNI_API __declspec(dllexport)
#elif defined(_WIN32)
#  define BNI_API __declspec(dllimport)
#elif defined(__GNUC__) || defined(__clang__)
#  define BNI_API __attribute__((visibility("default")))
#else
#  define BNI_API
#endif

#define BNI_VERSION_STRING "0.0.1"

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

typedef struct {
    int threads;
    int no_header_check;
    int force;
} bni_build_options_t;

typedef struct {
    uint64_t n_names;
    uint64_t n_records;
    uint64_t strings_size;
} bni_build_stats_t;

typedef struct {
    int threads;
    int ignore_metadata;
} bni_reader_options_t;

typedef struct bni_reader_t bni_reader_t;
typedef int (*bni_record_callback)(const bam1_t *record, const sam_hdr_t *header, void *user);

BNI_API int bni_write_index_file(const char *path, const bni_file_header_t *header,
                                 const bni_entry_t *entries, const char *strings);
BNI_API int bni_load_index_file(const char *path, bni_index_t *idx);
BNI_API void bni_index_destroy(bni_index_t *idx);
BNI_API const char *bni_entry_name(const bni_index_t *idx, const bni_entry_t *entry);
BNI_API const bni_entry_t *bni_find_entry(const bni_index_t *idx, const char *name);

BNI_API bni_index_t *bni_index_open(const char *path);
BNI_API void bni_index_close(bni_index_t *idx);

BNI_API int bni_build_index(const char *bam_path, const char *index_path,
                            const bni_build_options_t *opts, bni_build_stats_t *stats);

BNI_API int bni_reader_open(const char *bam_path, const char *index_path,
                            const bni_reader_options_t *opts, bni_reader_t **reader_out);
BNI_API void bni_reader_close(bni_reader_t *reader);
BNI_API const bni_index_t *bni_reader_index(const bni_reader_t *reader);
BNI_API const sam_hdr_t *bni_reader_header(const bni_reader_t *reader);
BNI_API int bni_reader_fetch(bni_reader_t *reader, const char *name,
                             bni_record_callback callback, void *user,
                             uint32_t *n_records_out);

#endif
