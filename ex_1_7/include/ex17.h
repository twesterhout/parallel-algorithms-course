#include <stdbool.h>
#include <stdint.h>

#if defined(_MSC_VER) || defined(WIN32) || defined(_WIN32)
#    define EX17_EXPORT __declspec(dllexport)
#    define EX17_NOINLINE __declspec(noinline)
#    define EX17_FORCEINLINE __forceinline inline
#    define EX17_LIKELY(cond) (cond)
#    define EX17_UNLIKELY(cond) (cond)
#    define EX17_CURRENT_FUNCTION __FUNCTION__
#else
#    define EX17_EXPORT __attribute__((visibility("default")))
#    define EX17_NOINLINE __attribute__((noinline))
#    define EX17_FORCEINLINE __attribute__((always_inline)) inline
#    define EX17_LIKELY(cond) __builtin_expect(!!(cond), 1)
#    define EX17_UNLIKELY(cond) __builtin_expect(!!(cond), 0)
#    define EX17_CURRENT_FUNCTION __PRETTY_FUNCTION__
#endif

typedef uint64_t ex17_prime_t;

typedef struct ex17_result {
    int           status;
    ex17_prime_t* primes;
    uint64_t      size;
} ex17_result_t;

ex17_result_t ex17_generate_primes_serial(ex17_prime_t upper_bound);
ex17_result_t ex17_generate_primes_parallel(ex17_prime_t upper_bound,
                                            uint64_t     block_size);

void ex17_filter_twins(uint64_t numbers[static 1], uint64_t* size);

int ex17_goldbach(uint64_t min, uint64_t max, bool* result);

#if defined(EX17_COUNTERS)
void     ex17_reset_scratch_counter(void);
uint64_t ex17_get_scratch_counter(void);
#endif
