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
    return getopt(argc, argv, "hqb:");
}

static void help(int argc, char* argv[])
{
    (void)argc;
    // clang-format off
    fprintf(stderr,
        "Usage: %s [-h] upper_bound\n"
        "\n"
        "Generates all prime numbers up to upper_bound and prints the number of\n"
        "'scratch-out' operations to stdout.\n", argv[0]);
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

    if (optind + 1 != argc) {
        fprintf(stderr, "%s: expected a single argument upper_bound\n",
                argv[0]);
        status = 1;
        goto cleanup;
    }
    uint64_t upper_bound;
    {
        char* end   = NULL;
        upper_bound = strtoull(argv[optind], &end, /*base=*/10);
        if (errno == ERANGE || (upper_bound == 0 && *end != '\0')) {
            fprintf(stderr, "%s: invalid upper_bound: '%s'\n", argv[0],
                    argv[optind]);
            status = 1;
            goto cleanup;
        }
    }

    ex17_reset_scratch_counter();
    ex17_result_t result;
    result = ex17_generate_primes_serial(upper_bound);

    if (result.status != 0) {
        fprintf(stderr, "ex17_generate_primes_serial() failed: %s\n",
                strerror(result.status));
        status = result.status;
        goto cleanup;
    }

    fprintf(stdout, "%lu\n", ex17_get_scratch_counter());
    if (result.primes != NULL) { free(result.primes); }

cleanup:
    return status;
}
