#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "./system_gen.h"

/*
  Computes a 2-norm of a vector v.
*/
static double vector_norm(double *v, int n) {

  double sum_square = 0.0;
  for (int i = 0; i < n; i++)
    sum_square += v[i] * v[i];

  return sqrt(sum_square);
}

/*
  Computes a 2-norm of Ax-b.
*/
static double residual_norm(double *A, double *x, double *b, int n) {

  double sum_square = 0.0;

  for (int i = 0; i < n; i++) {
    double Ax_i = 0.0;

    for (int j = 0; j < n; j++)
      Ax_i += A[i*n + j] * x[j];

    double r = Ax_i - b[i];
    sum_square += r*r;
  }

  return sqrt(sum_square);
}

/*
  Solves Ax=b via the Jacobi method with an initial solution in x.

  Inputs
    Max_iter: maximum number of Jacobi iterations to try.
    tolerance: relative residual error used to declare convergence.
    x: the initial solution vector of size n.

  Output:
    x: the solution to Ax=b.
*/
static double *jacobi(double *A, double *b, double *x,
                      int n, int max_iter, double tolerance) {

  double *x_new = calloc(n, sizeof(double));
  double b_norm = vector_norm(b, n);

  int iter;
  for (iter = 0; iter < max_iter; iter++) {
    for (int i = 0; i < n; i++) {
      double sigma = 0.0;
      for (int j = 0; j < n; j++)
        if (j != i)
          sigma += A[i*n + j] * x[j];
      x_new[i] = (b[i] - sigma) / A[i*n + i];
    }
    double *tmp = x;
    x = x_new;
    x_new = tmp;

    // If the solution converged, stop iterating.
    double accuracy = residual_norm(A, x, b, n);
    if (accuracy / b_norm < tolerance) {
      printf("Solved with error %e in %d iterations. \n",
             accuracy / b_norm, iter);
      break;
    }
  }

  free(x_new);
  if (iter == max_iter)
    printf("Solver did not converge in %d iterations. \n",
           max_iter);
  return x;
}

int main(int argc, char *argv[]) {

  if (argc < 3) {
    fprintf(stderr,
      "Please specify n the size of the square linear system and margin.\n");
    return 1;
  }

  int n = atoi(argv[1]);
  double margin = atof(argv[2]);

  double *A = calloc((size_t) n * n, sizeof(double));
  double *b = calloc(n, sizeof(double));
  double *x = calloc(n, sizeof(double));

  if (!A || !b || !x) {
    fprintf(stderr, "Memory allocation for A, b, or x has failed. \n");
    return 1;
  }

  double t0 = clock();

  // Generate a linear system with diagonally dominant matrix A.
  generate_system(A, b, n, 0, margin);

  double tolerance = 1e-10;
  int max_iter = 500;

  // Attempt to solve the linear system.
  x = jacobi(A, b, x, n, max_iter, tolerance);

  double runTime = (clock() - t0)/(double)CLOCKS_PER_SEC;
  printf("Total run time:     %f seconds\n", runTime);

  free(A);
  free(b);
  free(x);
}
