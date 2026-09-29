#ifndef FFT_H
#define FFT_H
#include <complex.h>
#include <stddef.h>

/* n must be a positive power of two. Both routines modify data in place. */
void fft1d(double complex *data, int n);
/* Included for reference. Not used in the assessed project.
   fft1d followed by ifft1d restores the input up to rounding error. */
void ifft1d(double complex *data, int n);

#define FFT_MAX_N 16384
#define FFT_REFERENCE_LIMIT 64
#define FFT_BENCHMARK_TOLERANCE 1e-10
int parse_positive_int(const char *text, int *value);
int valid_fft_size(int n);
void *fft_allocate(size_t count);
void generate_input(double complex *data, int n);
void transpose_serial(const double complex *input, double complex *output, int n);
void fft2d_serial(double complex *data, double complex *scratch, int n);
double max_error(const double complex *actual, const double complex *expected,
                 size_t count);
double benchmark_error(const double complex *actual, int n, int first_row,
                       int rows);
void print_matrix(const double complex *data, int n);
#endif
