#include "ex17.h"

#include <getopt.h>
#include <unistd.h>

#include <ctype.h> // isprint
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static inline int _parse(int argc, char* argv[])
{
    return getopt(argc, argv, "h");
}

static void help(int argc, char* argv[])
{
    (void)argc;
    // clang-format off
    fprintf(stderr,
        "Usage: %s [-h] min max\n"
        "\n"
        "Checks the Goldbach conjecture for all integers in the range [min, max).\n",
        argv[0]);
    // clang-format on
}

int main(int argc, char* argv[])
{
    opterr         = 0;
    int  status    = 0;
    bool help_flag = false;
    for (int c = _parse(argc, argv); c != -1; c = _parse(argc, argv)) {
        switch (c) {
        case 'h': help_flag = true; break;
        case '?':
            if (isprint(optopt)) {
                fprintf(stderr,
                        "Unknown option '-%c'. Use -h flag to see "
                        "available options\n",
                        optopt);
            }
            else {
                fprintf(stderr,
                        "Unknown option character '\\x%x'. Use -h flag to "
                        "see available options\n",
                        optopt);
            }
            status = 1;
            goto cleanup;
        } // end switch
    }
    if (help_flag) {
        help(argc, argv);
        goto cleanup;
    }

    if (optind + 2 != argc) {
        fprintf(stderr, "%s: expected two arguments: min and max\n", argv[0]);
        status = 1;
        goto cleanup;
    }
    uint64_t min;
    {
        char* end = NULL;
        min       = strtoull(argv[optind], &end, /*base=*/10);
        if (errno == ERANGE || (min == 0 && *end != '\0') || min < 3) {
            fprintf(stderr, "%s: invalid min: '%s'\n", argv[0], argv[optind]);
            status = 1;
            goto cleanup;
        }
    }
    uint64_t max;
    {
        char* end = NULL;
        max       = strtoull(argv[optind + 1], &end, /*base=*/10);
        if (errno == ERANGE || (max == 0 && *end != '\0') || max < min) {
            fprintf(stderr, "%s: invalid max: '%s'\n", argv[0],
                    argv[optind + 1]);
            status = 1;
            goto cleanup;
        }
    }

    bool result;
    status = ex17_goldbach(min, max, &result);

    if (status != 0) {
        fprintf(stderr, "ex17_goldbach() failed: %s\n", strerror(status));
        goto cleanup;
    }

    printf("%s\n", result ? "true" : "false");

cleanup:
    return status;
}
