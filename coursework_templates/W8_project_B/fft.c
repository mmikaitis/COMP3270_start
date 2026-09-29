#include "fft.h"
#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/* Supplied numerical kernel. Students do not need to modify this file. */
static void transform(double complex *data, int n, int inverse) {
    const double pi = 3.14159265358979323846;
    assert(n > 0 && (n & (n - 1)) == 0);
    /* Put the input into bit-reversed order. */
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            double complex value = data[i];
            data[i] = data[j];
            data[j] = value;
        }
    }
    /* Combine transforms of length 2, 4, 8, ... */
    for (int length = 2; length <= n; length *= 2) {
        double angle = (inverse ? 2.0 : -2.0) * pi / length;
        double complex step = cos(angle) + I * sin(angle);
        for (int start = 0; start < n; start += length) {
            double complex factor = 1.0;
            for (int j = 0; j < length / 2; ++j) {
                double complex a = data[start + j];
                double complex b = factor * data[start + j + length / 2];
                data[start + j] = a + b;
                data[start + j + length / 2] = a - b;
                factor *= step;
            }
        }
        if (length == n) break; /* Avoid overflow in the loop increment. */
    }
    if (inverse)
        for (int i = 0; i < n; ++i) data[i] /= n;
}
void fft1d(double complex *data, int n) { transform(data, n, 0); }
void ifft1d(double complex *data, int n) { transform(data, n, 1); }


int parse_positive_int(const char *text, int *value) {
    char *end;
    errno = 0;
    long parsed = strtol(text, &end, 10);
    if (errno || end == text || *end || parsed < 1 || parsed > INT_MAX) return 0;
    *value = (int)parsed;
    return 1;
}
int valid_fft_size(int n) {
    return n >= 1 && n <= FFT_MAX_N && (n & (n - 1)) == 0;
}
void *fft_allocate(size_t count) {
    if (!count || count > SIZE_MAX / sizeof(double complex)) return NULL;
    return calloc(count, sizeof(double complex));
}
void generate_input(double complex *data, int n) {
    const double pi = 3.14159265358979323846;
    uint32_t state = 1234567;
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            if (n <= FFT_REFERENCE_LIMIT) {
                /* General, nonsymmetric complex values, generated once on rank 0. */
                state = 1664525u * state + 1013904223u;
                double re = (double)state / 4294967296.0 - 0.5;
                state = 1664525u * state + 1013904223u;
                double im = (double)state / 4294967296.0 - 0.5;
                data[(size_t)r*n+c] = re + I*im;
            } else {
                /* Supplied benchmark data with a known spectrum.
                   Students do not need to derive these expressions. */
                double a = 2*pi*(2*r+3*c)/n, b = 2*pi*(5*r+c)/n;
                data[(size_t)r*n+c] =
                    cos(a) + I*sin(a) + 0.5*(cos(b) + I*sin(b));
            }
        }
    }
}
void transpose_serial(const double complex *input, double complex *output, int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            output[(size_t)j*n+i] = input[(size_t)i*n+j];
}
void fft2d_serial(double complex *data, double complex *scratch, int n) {
    /* Same numerical routine and per-row arithmetic as the MPI version. */
    for (int row = 0; row < n; ++row) fft1d(data + (size_t)row*n, n);
    transpose_serial(data, scratch, n);
    for (int row = 0; row < n; ++row) fft1d(scratch + (size_t)row*n, n);
    transpose_serial(scratch, data, n);
}
double max_error(const double complex *actual, const double complex *expected,
                 size_t count) {
    double largest = 0;
    for (size_t i = 0; i < count; ++i) {
        double error = cabs(actual[i] - expected[i]);
        if (!isfinite(error)) return INFINITY;
        if (error > largest) largest = error;
    }
    return largest;
}
double benchmark_error(const double complex *actual, int n, int first_row,
                       int rows) {
    double largest = 0, scale = (double)n*n;
    for (int i = 0; i < rows; ++i) {
        int r = first_row + i;
        for (int c = 0; c < n; ++c) {
            double expected = 0;
            if (r == 2 && c == 3) expected = scale;
            if (r == 5 && c == 1) expected = 0.5*scale;
            double error = cabs(actual[(size_t)i*n+c] - expected) / scale;
            if (!isfinite(error)) return INFINITY;
            if (error > largest) largest = error;
        }
    }
    return largest;
}
void print_matrix(const double complex *data, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double complex value = data[(size_t)i*n+j];
            printf(" % .6g%+.6gi", creal(value), cimag(value));
        }
        putchar('\n');
    }
}
