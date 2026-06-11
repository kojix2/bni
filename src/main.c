#define _POSIX_C_SOURCE 200809L

#include "bni.h"

#include <stdio.h>
#include <string.h>

static void usage(FILE *fp) {
    fprintf(fp,
        "bni %s - BAM Name Index for queryname-sorted BAM files\n\n"
        "Usage:\n"
        "  bni <command> [options]\n\n"
        "Commands:\n"
        "  index     build a .bam.bni index for a samtools sort -N BAM\n"
        "  get       extract records by read name\n"
        "  stats     show BNI index metadata\n"
        "  check     verify BAM/index consistency\n"
        "  version   print version\n\n"
        "Run 'bni <command> --help' for command-specific help.\n",
        BNI_VERSION_STRING);
}

int main(int argc, char **argv) {
    if (argc < 2) { usage(stderr); return 1; }
    if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) { usage(stdout); return 0; }
    if (strcmp(argv[1], "version") == 0 || strcmp(argv[1], "--version") == 0) { printf("bni %s\n", BNI_VERSION_STRING); return 0; }
    if (strcmp(argv[1], "index") == 0) return bni_cmd_index(argc - 1, argv + 1);
    if (strcmp(argv[1], "get") == 0 || strcmp(argv[1], "query") == 0) return bni_cmd_get(argc - 1, argv + 1);
    if (strcmp(argv[1], "stats") == 0) return bni_cmd_stats(argc - 1, argv + 1);
    if (strcmp(argv[1], "check") == 0) return bni_cmd_check(argc - 1, argv + 1);
    bni_print_error("unknown command '%s'", argv[1]);
    usage(stderr);
    return 1;
}
