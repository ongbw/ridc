#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <stdio.h>

#include "ridc.h"

/** @brief This is the main function for the implicit_mkl example
 *
 * This will pass user given options along with some standard options
 * for this type of problem into the PARAMITER struct and start 
 * the solving process by calling ridc_be()
 * \n
 * This example uses the Math Kernal Library (mlk) 
 */

int main(int argc, char *argv[]) {
  int order, nt;
  double *sol;

  if (argc != 3) {
    printf("usage: <executable> <order> <nt>  >  output_file\n");
    fflush(stdout);
    exit(1);
  }
  else {
    order = atoi(argv[1]); // order of method
    nt = atoi(argv[2]); // number of time steps
  }


  
  // initialize parameters for problem
  PARAMETER param;
  param.ti = 0; // initial time
  param.tf = 1; // final time
  param.dt = (double)(param.tf - param.ti)/nt; // compute dt
  param.neq = 2; // number of equations
  param.nt = nt; // store number of time steps


  sol = new double[param.neq];
  // specify initial condition
  for (int i =0; i<param.neq; i++) {
    sol[i]=1.0;
  }


  // call ridc 
  ridc_be(order, param, sol);

  // output solution to screen
  for (int i = 0; i < param.neq; i++)
    printf("%14.12f\n", sol[i]);
  delete [] sol;
}
