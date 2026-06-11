#ifndef BNI_INTERNAL_H
#define BNI_INTERNAL_H

#include "bni.h"

#include <stddef.h>
#include <stdint.h>

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
