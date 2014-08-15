#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <stdio.h>

#include "ridc.h"
#include "ode.h"

int main(int argc, char *argv[]) {
  int order, nt;
  double *sol;

  if (argc != 3) {
    printf("usage: <executable> <order> <nt>  >  output_file\n");
    fflush(stdout);
    exit(1);
  } else {
    order = atoi(argv[1]); // order of method
    nt = atoi(argv[2]); // number of time steps
  }

  // initialize parameters from problem
  PARAMETER param;
  init_params(nt, param);

  sol = new double[param.neq];

  // select appropriate functions based on compilation flag
#if EXPLICIT
  ridc_fe(order, param, sol);
#elif IMPLICIT
  ridc_be(order, param, sol);
#else
#error Neither IMPLICIT nor EXPLICIT flags are set!
#endif

  // output solution to screen
  for (int i = 0; i < param.neq; i++)
    printf("%14.12f\n", sol[i]);
  delete [] sol;
}
