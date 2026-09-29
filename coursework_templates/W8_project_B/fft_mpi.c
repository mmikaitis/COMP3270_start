#include "fft.h"
#include <mpi.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* One rank owns m consecutive rows, each containing n complex values.
   All ranks follow the same sequence of collective calls. */
static void fail(const char *message) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) fprintf(stderr, "%s\n", message);
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
}

static double complex *allocate_or_abort(size_t count) {
    double complex *data = fft_allocate(count);
    if (!data) fail("Array allocation failed");
    return data;
}

/* Used for the initial input and for small reference arrays. */
static void scatter_rows(const double complex *full, double complex *local,
                         int count) {
    /* TODO Task 2: distribute count complex values to each rank.
       full is valid only on rank 0. Complete the MPI_Scatter call. */
    (void)full; (void)local; (void)count;
    fail("TODO Task 2: complete scatter_rows");
}

static void distributed_transpose(const double complex *input,
                                  double complex *output,
                                  double complex *send, double complex *recv,
                                  int n, int m, int processes) {
    size_t block_count = (size_t)m*m;
    (void)block_count; (void)n; /* Used when you complete the indices and exchange. */
    /* Each destination receives one m x m block.
       Blocks are stored by destination, then by row within each block. */
    for (int q = 0; q < processes; ++q) {
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j) {
                /* TODO Task 2: fill both indices, then remove this fail call.
                   input has m rows of length n. Destination q needs its m columns. */
                fail("TODO Task 2: complete the packing indices");
                size_t source = 0;
                size_t packed = 0;
                send[packed] = input[source];
            }
        }
    }

    /* TODO Task 2: exchange one block per destination with MPI_Alltoall.
       Use block_count as the count per destination, not the total buffer size. */
    fail("TODO Task 2: complete the all-to-all exchange");

    /* Received blocks are ordered by source rank.
       Place each value into a complete row of the transposed matrix. */
    for (int r = 0; r < processes; ++r) {
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j) {
                /* TODO Task 2: fill both indices, then remove this fail call.
                   Block r came from source rank r. Transpose within the block. */
                fail("TODO Task 2: complete the unpacking indices");
                size_t received = 0;
                size_t destination = 0;
                output[destination] = recv[received];
            }
        }
    }
}

static double global_error(double local_error) {
    double error = INFINITY;
    /* TODO Task 2: use MPI_Allreduce to find the largest local_error.
       Each rank must receive the result. The scalar datatype is MPI_DOUBLE. */
    (void)local_error;
    fail("TODO Task 2: complete the global error reduction");
    return error;
}

static double slowest_time(double local_time) {
    double elapsed = INFINITY;
    /* TODO Task 3: use MPI_Allreduce to report the slowest rank's time. */
    (void)local_time;
    fail("TODO Task 3: complete the timing reduction");
    return elapsed;
}

/* Supplied stage sequence. The result returns to ordinary row layout. */
static void fft2d_distributed(double complex *data, double complex *scratch,
                              double complex *send, double complex *recv,
                              int n, int m, int processes) {
    for (int row = 0; row < m; ++row) fft1d(data + (size_t)row*n, n);
    distributed_transpose(data, scratch, send, recv, n, m, processes);
    for (int row = 0; row < m; ++row) fft1d(scratch + (size_t)row*n, n);
    distributed_transpose(scratch, data, send, recv, n, m, processes);
}

static void check_transpose(double complex *data, double complex *scratch,
                            double complex *send, double complex *recv,
                            int n, int m, int rank, int processes) {
    /* Generate owned rows of A[i,j] = 10*i+j. No scatter needed for this test. */
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            data[(size_t)i*n+j] = 10.0*(rank*m+i) + j;
    distributed_transpose(data, scratch, send, recv, n, m, processes);
    double local_error = 0;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            double expected = 10.0*j + rank*m+i;
            double error = cabs(scratch[(size_t)i*n+j] - expected);
            if (!isfinite(error)) local_error = INFINITY;
            else if (error > local_error) local_error = error;
        }
    }
    double first_error = global_error(local_error);
    distributed_transpose(scratch, data, send, recv, n, m, processes);
    local_error = 0;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            double error = cabs(data[(size_t)i*n+j] - (10.0*(rank*m+i)+j));
            if (!isfinite(error)) local_error = INFINITY;
            else if (error > local_error) local_error = error;
        }
    }
    double restored_error = global_error(local_error);
    if (rank == 0)
        printf("Transpose N=%d ranks=%d: error=%.17g twice=%.17g\n",
               n, processes, first_error, restored_error);
    if (first_error != 0 || restored_error != 0) fail("Transpose check failed");
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    int n = 0, repetitions = 0;
    int transpose_only = argc >= 2 && strcmp(argv[1], "--transpose") == 0;
    int check_only = argc >= 2 && strcmp(argv[1], "--check") == 0;
    int valid = argc == 3;
    if (valid && (transpose_only || check_only))
        valid = parse_positive_int(argv[2], &n);
    else if (valid)
        valid = parse_positive_int(argv[1], &n) &&
                parse_positive_int(argv[2], &repetitions) && repetitions <= 100;
    valid = valid && valid_fft_size(n) && processes <= n && n % processes == 0;
    if (!valid) {
        if (rank == 0)
            fprintf(stderr, "Usage: %s N repetitions | --transpose N | --check N\n"
                    "N: power of two, 1..%d, divisible by the rank count.\n"
                    "Repetitions: 1..100. Aire benchmark: N=16384, 8..64 ranks, 3 repeats.\n",
                    argv[0], FFT_MAX_N);
        MPI_Finalize();
        return 1;
    }
    int m = n / processes;
    size_t count = (size_t)m*n;
    if (count > INT_MAX || (size_t)m*m > INT_MAX)
        fail("Array is too large for MPI counts");
    double complex *data = allocate_or_abort(count);
    double complex *scratch = allocate_or_abort(count);
    double complex *send = allocate_or_abort(count);
    double complex *recv = allocate_or_abort(count);
    if (transpose_only) {
        check_transpose(data, scratch, send, recv, n, m, rank, processes);
        free(data); free(scratch); free(send); free(recv);
        MPI_Finalize();
        return 0;
    }

    double complex *original = allocate_or_abort(count);
    double complex *full = NULL;
    if (rank == 0) {
        full = allocate_or_abort((size_t)n*n);
        generate_input(full, n);
    }
    scatter_rows(full, original, (int)count);

    /* Small tests compare with the same sequential FFT.
       Larger runs use known benchmark coefficients instead of a root FFT. */
    double complex *reference = NULL;
    if (n <= FFT_REFERENCE_LIMIT) {
        reference = allocate_or_abort(count);
        if (rank == 0) {
            double complex *temporary = allocate_or_abort((size_t)n*n);
            fft2d_serial(full, temporary, n);
            free(temporary);
        }
        scatter_rows(full, reference, (int)count);
    }
    free(full);
    double tolerance = reference ? 0.0 : FFT_BENCHMARK_TOLERANCE;

    /* repeat=0 is the warm-up. All ranks reset input before each transform. */
    for (int repeat = 0; repeat <= repetitions; ++repeat) {
        memcpy(data, original, count*sizeof(*data));
        if (!check_only) MPI_Barrier(MPI_COMM_WORLD);
        double start = check_only ? 0 : MPI_Wtime();
        fft2d_distributed(data, scratch, send, recv, n, m, processes);
        double local_time = check_only ? 0 : MPI_Wtime() - start;
        double elapsed = check_only ? 0 : slowest_time(local_time);
        double local_error = reference ? max_error(data, reference, count) :
                                        benchmark_error(data, n, rank*m, m);
        double error = global_error(local_error);
        if (!isfinite(error) || error > tolerance) {
            if (rank == 0) {
                fprintf(stderr,
                        "FFT check failed: N=%d ranks=%d repeat=%d "
                        "error=%.17g tolerance=%.17g\n",
                        n, processes, repeat, error, tolerance);
                fflush(stderr);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            fail("FFT check failed");
        }
        if (!check_only && (!isfinite(elapsed) || elapsed <= 0))
            fail("Invalid elapsed time");
        if (check_only && rank == 0)
            printf("Check N=%d ranks=%d: error=%.17g\n", n, processes, error);
        if (repeat > 0 && rank == 0)
            printf("%d %d %.9g %.17g\n", processes, repeat, elapsed, error);
    }

    if (rank == 0)
        fprintf(stderr, "PASS N=%d ranks=%d (%s reference)\n",
                n, processes, reference ? "exact serial" : "scaled mathematical");
    free(reference); free(original);
    free(data); free(scratch); free(send); free(recv);
    MPI_Finalize();
    return 0;
}
