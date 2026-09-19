#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Uniformly distributed random numbers in the range (0, 1). */
static double rand_uniform_open01(void) {
  return ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
}

/* Normal distribution random numbers */
static double rand_normal_standard(void) {
  /* Box-Muller transform */
  const double PI = 3.14159265358979323846;
  double u1 = rand_uniform_open01();
  double u2 = rand_uniform_open01();
  return sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
}

void simulateQueue(double T, double lambda_arrival, double mu_service, double C,
                   int runs) {

  const double dt = 0.01;   /* timestep [0.01 minutes = 0.6 seconds] */
  const int N = (int)(T / dt);
  const double sigma = 1.0; /* noise strength (stdev of 1 job/min) */

  double Qsum = 0.0;
  double Qdev = 0.0;

  for (int r = 0; r < runs; r++) {

    /* Allocate and initialise storage for Q, N timesteps. */
    double *Q = (double *)calloc((size_t)N, sizeof(double));
    if (Q == NULL) {
      fprintf(stderr, "Memory allocation failed.\n");
      return;
    }

    /* Update the model for N timesteps. */
    for (int t = 0; t < N - 1; t++) {
      double q_min_c = (Q[t] < C) ? Q[t] : C;
      double drift = lambda_arrival - mu_service * q_min_c;
      double diffusion = sigma * sqrt(dt) * rand_normal_standard();

      Q[t + 1] = Q[t] + drift * dt + diffusion;
      if (Q[t + 1] < 0.0) {
        Q[t + 1] = 0.0;
      }
    }

    double mean = 0.0;
    for (int i = 0; i < N; i++) {
      mean += Q[i];
    }
    mean /= (double)N;

    /* Population variance. */
    double var = 0.0;
    for (int i = 0; i < N; i++) {
      double d = Q[i] - mean;
      var += d * d;
    }
    var /= (double)N;

    /* Compute sums of means and standard deviations. */
    Qsum += mean;
    Qdev += sqrt(var);

    free(Q);
  }

  printf("Queue statistics over %d runs\n", runs);
  printf("-------------------------------\n");
  printf("Mean queue length:  %f\n", Qsum / (double)runs);
  printf("Stdev queue length: %f\n", Qdev / (double)runs);
  printf("Mean occupancy:     %f%%\n", 100.0 * Qsum / ((double)runs * C));
}

int main(int argc, char *argv[]) {
  int runs = 1;

  if (argc > 1) {
    runs = atoi(argv[1]);
    if (runs <= 0) {
      fprintf(stderr, "Invalid runs value. Using default runs = 1.\n");
      runs = 1;
    }
  }

  /* Set up random number generator seed */
  srand((unsigned int)time(NULL));

  simulateQueue(120.0, 400.0 / 120.0, 0.2, 32.0, runs);

  return 0;
}
