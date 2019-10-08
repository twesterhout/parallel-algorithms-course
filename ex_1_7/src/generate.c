#include "ex17.h"

#include <mpi.h>

#include <getopt.h>
#include <unistd.h>

#include <assert.h>
#include <ctype.h> // isprint
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static inline int _parse(int argc, char* argv[])
{
    return getopt(argc, argv, "hqtb:");
}

static void help(int argc, char* argv[])
{
    (void)argc;
    // clang-format off
    fprintf(stderr,
        "Usage: %s [-q] [-t] [-b block_size] upper_bound\n"
        "\n"
        "Generates all prime numbers up to upper_bound and prints them to stdout.\n"
        "If block_size is given, parallel algorithm is used, otherwise primes are\n"
        "generated using a serial implementation.\n"
        "\n"
        "-q flag causes output to be supressed. This is useful for benchmarking.\n"
        "-t flag causes the program to generate twin prime numbers.\n", argv[0]);
    // clang-format on
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int status = 0;

    // Make sure getopt doesn't try to print error messages from all the
    // processes
    opterr               = 0;
    bool  help_flag      = false;
    bool  quiet_flag     = false;
    bool  twins_flag     = false;
    char* block_size_str = NULL;
    for (int c = _parse(argc, argv); c != -1; c = _parse(argc, argv)) {
        switch (c) {
        case 'b': block_size_str = optarg; break;
        case 'q': quiet_flag = true; break;
        case 't': twins_flag = true; break;
        case 'h': help_flag = true; break;
        case '?':
            if (rank == 0) {
                if (optopt == 'b') {
                    fprintf(stderr, "Option -%c requires an argument.\n",
                            optopt);
                }
                else if (isprint(optopt)) {
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
            }
            status = 1;
            goto cleanup;
        } // end switch
    }
    if (help_flag) {
        if (rank == 0) { help(argc, argv); }
        goto cleanup;
    }

    if (optind + 1 != argc) {
        if (rank == 0) {
            fprintf(stderr, "%s: expected a single argument upper_bound\n",
                    argv[0]);
        }
        status = 1;
        goto cleanup;
    }
    uint64_t upper_bound;
    {
        char* end   = NULL;
        upper_bound = strtoull(argv[optind], &end, /*base=*/10);
        if (errno == ERANGE || (upper_bound == 0 && *end != '\0')) {
            if (rank == 0) {
                fprintf(stderr, "%s: invalid upper_bound: '%s'\n", argv[0],
                        argv[optind]);
            }
            status = 1;
            goto cleanup;
        }
    }

    ex17_result_t result;
    if (block_size_str != NULL) {
        char*    end        = NULL;
        uint64_t block_size = strtoull(block_size_str, &end, /*base=*/10);
        if (errno == ERANGE || block_size == 0 || *end != '\0') {
            if (rank == 0) {
                fprintf(stderr, "%s: invalid block_size: '%s'\n", argv[0],
                        block_size_str);
            }
            status = 1;
            goto cleanup;
        }
        result = ex17_generate_primes_parallel(upper_bound, block_size);
    }
    else {
        result = ex17_generate_primes_serial(upper_bound);
    }

    if (result.status != 0) {
        if (rank == 0) {
            fprintf(stderr, "ex17_generate_primes() failed: %s\n",
                    strerror(result.status));
        }
        status = result.status;
        MPI_Abort(MPI_COMM_WORLD, status);
        goto cleanup;
    }
    if (!quiet_flag && rank == 0) {
        if (twins_flag) {
            ex17_filter_twins(result.primes, &result.size);
            assert(result.size % 2 == 0);
            for (uint64_t i = 0; i < result.size; i += 2) {
                printf("%lu %lu\n", result.primes[i], result.primes[i + 1]);
            }
        }
        else {
            for (uint64_t i = 0; i < result.size; ++i) {
                printf("%lu\n", result.primes[i]);
            }
        }
    }
    if (result.primes != NULL) { free(result.primes); }

cleanup:
    MPI_Finalize();
    return status;
}
