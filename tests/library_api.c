#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <bni.h>

typedef struct {
    const char *expected;
    uint32_t count;
} count_context_t;

static int count_record(const bam1_t *record, const sam_hdr_t *header, void *user) {
    (void)header;
    count_context_t *ctx = (count_context_t *)user;
    const char *qname = bam_get_qname(record);
    if (qname == NULL || strcmp(qname, ctx->expected) != 0) return -1;
    ctx->count++;
    return 0;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "usage: library_api <bam> <qname> <expected-count>\n");
        return 2;
    }
    uint32_t expected_count = (uint32_t)strtoul(argv[3], NULL, 10);
    bni_reader_t *reader = NULL;
    if (bni_reader_open(argv[1], NULL, NULL, &reader) != 0) return 1;

    count_context_t ctx;
    ctx.expected = argv[2];
    ctx.count = 0;
    uint32_t fetched_count = 0;
    int rc = bni_reader_fetch(reader, argv[2], count_record, &ctx, &fetched_count);
    bni_reader_close(reader);

    if (rc != 0) return 1;
    if (ctx.count != expected_count || fetched_count != expected_count) {
        fprintf(stderr, "expected %u records, callback=%u fetched=%u\n",
                expected_count, ctx.count, fetched_count);
        return 1;
    }
    return 0;
}
