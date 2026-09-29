#include "fft.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int run_example(const char *name, const double complex *input,
                       const double complex *expected) {
    double complex data[4];
    for (int i = 0; i < 4; ++i) data[i] = input[i];
    /* TODO Task 1: call the supplied forward FFT on data (length 4).
       Remove the next two lines when you add your call. */
    fprintf(stderr, "TODO Task 1: call fft1d in run_example\n");
    return 1;
    printf("%s:", name);
    for (int i = 0; i < 4; ++i)
        printf(" %g%+gi", creal(data[i]), cimag(data[i]));
    double error = max_error(data, expected, 4);
    printf("\nError: %.17g\n", error);
    return error != 0;
}
static int run_examples(void) {
    const double complex constant[4] = {1,1,1,1};
    const double complex constant_fft[4] = {4,0,0,0};
    const double complex alternating[4] = {1,-1,1,-1};
    const double complex alternating_fft[4] = {0,0,4,0};
    int failed = run_example("Constant", constant, constant_fft);
    failed |= run_example("Alternating", alternating, alternating_fft);
    if (failed) fprintf(stderr, "Example check failed\n");
    return failed;
}

int main(int argc, char **argv) {
    if (argc == 2 && strcmp(argv[1], "--examples") == 0)
        return run_examples();
    int n;
    if (argc != 2 || !parse_positive_int(argv[1], &n) || !valid_fft_size(n)) {
        fprintf(stderr, "Usage: %s N (power of two, 1..%d) or --examples\n", argv[0], FFT_MAX_N);
        return 1;
    }
    size_t count = (size_t)n*n;
    double complex *data = fft_allocate(count), *scratch = fft_allocate(count);
    if (!data || !scratch) {
        fprintf(stderr, "Allocation failed\n");
        free(data); free(scratch);
        return 1;
    }
    generate_input(data, n);
    if (n <= 16) { puts("Input:"); print_matrix(data, n); }
    fft2d_serial(data, scratch, n);
    if (n <= 16) { puts("FFT output:"); print_matrix(data, n); }
    int failed = 0;
    if (n > FFT_REFERENCE_LIMIT) {
        double error = benchmark_error(data, n, 0, n);
        printf("Scaled error against known spectrum: %.17g\n", error);
        failed = error > FFT_BENCHMARK_TOLERANCE;
    }
    printf("Completed serial %d x %d FFT\n", n, n);
    free(data); free(scratch);
    return failed;
}
