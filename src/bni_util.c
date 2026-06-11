#define _POSIX_C_SOURCE 200809L

#include "bni.h"

#include <errno.h>
#include <limits.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

char *bni_strdup(const char *s) {
    if (s == NULL) return NULL;
    size_t n = strlen(s) + 1;
    char *p = (char *)malloc(n);
    if (p == NULL) return NULL;
    memcpy(p, s, n);
    return p;
}

char *bni_default_index_path(const char *bam_path) {
    if (bam_path == NULL) return NULL;
    size_t n = strlen(bam_path);
    const char *suffix = ".bni";
    size_t m = strlen(suffix);
    if (n > SIZE_MAX - m - 1) return NULL;
    char *out = (char *)malloc(n + m + 1);
    if (out == NULL) return NULL;
    memcpy(out, bam_path, n);
    memcpy(out + n, suffix, m + 1);
    return out;
}

uint64_t bni_fnv1a64(const void *data, size_t len) {
    const unsigned char *p = (const unsigned char *)data;
    uint64_t h = UINT64_C(14695981039346656037);
    for (size_t i = 0; i < len; ++i) {
        h ^= (uint64_t)p[i];
        h *= UINT64_C(1099511628211);
    }
    return h;
}

int bni_file_metadata(const char *path, uint64_t *size_out, int64_t *mtime_out) {
    struct stat st;
    if (stat(path, &st) != 0) return -1;
    if (size_out) *size_out = (uint64_t)st.st_size;
    if (mtime_out) *mtime_out = (int64_t)st.st_mtime;
    return 0;
}

int bni_path_exists(const char *path) {
    return path != NULL && access(path, F_OK) == 0;
}

int bni_has_suffix(const char *s, const char *suffix) {
    if (s == NULL || suffix == NULL) return 0;
    size_t n = strlen(s), m = strlen(suffix);
    return m <= n && strcmp(s + n - m, suffix) == 0;
}

int bni_parse_threads(const char *s, int *out) {
    if (s == NULL || *s == '\0') return -1;
    errno = 0;
    char *end = NULL;
    long v = strtol(s, &end, 10);
    if (errno != 0 || *end != '\0' || v < 0 || v > INT_MAX) return -1;
    if (out) *out = (int)v;
    return 0;
}

void bni_print_error(const char *fmt, ...) {
    va_list ap;
    fprintf(stderr, "bni: error: ");
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

void bni_print_warning(const char *fmt, ...) {
    va_list ap;
    fprintf(stderr, "bni: warning: ");
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

void bni_format_u64(char *buf, size_t buflen, uint64_t value) {
    char tmp[64];
    int pos = 0, group = 0;
    if (buflen == 0) return;
    if (value == 0) { snprintf(buf, buflen, "0"); return; }
    while (value > 0 && pos < (int)sizeof(tmp) - 1) {
        if (group == 3) { tmp[pos++] = ','; group = 0; }
        tmp[pos++] = (char)('0' + (value % 10));
        value /= 10;
        group++;
    }
    size_t out = 0;
    while (pos > 0 && out + 1 < buflen) buf[out++] = tmp[--pos];
    buf[out] = '\0';
}
