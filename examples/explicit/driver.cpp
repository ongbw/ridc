#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <stdio.h>

#include "ridc.h"

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


  
  // initialize parameters for problem
  PARAMETER param;
  param.ti = 0;
  param.tf = 1;
  param.dt = (double)(param.tf - param.ti)/nt;
  param.neq = 2; // some neq
  param.nt = nt;


  sol = new double[param.neq];
  // initial condition
  for (int i =0; i<param.neq; i++) {
    sol[i]=1.0;
  }


  ridc_fe(order, param, sol);

  // output solution to screen
  for (int i = 0; i < param.neq; i++)
    printf("%14.12f\n", sol[i]);
  delete [] sol;
}
