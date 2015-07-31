#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <stdio.h>
#include "mkl.h" // MKL Standard Lib
#include "mkl_lapacke.h" // MKL LAPACK
#include <cmath>

#include "ridc.h"

//EditNoticeForOngbw added brief and detailed descripion below for Doxygen to use

/** @brief This is the main function for the brusselator_mkl example
 *
 * This will pass user given options along with some standard options
 * for this type of problem into the PARAMITER struct and start 
 * the solving process by calling ridc_be()
 * \n
 * This example uses the Math Kernal Library (mlk) 
 */
 
int main(int argc, char *argv[]) {
  int order, nt,neq;
  double *sol;

  if (argc != 4) {
    printf("usage: <executable> <order> <nt>  >  output_file\n");
    fflush(stdout);
    exit(1);
  }
  else {
    order = atoi(argv[1]); // order of method
    nt = atoi(argv[2]); // number of time steps
    neq = atoi(argv[3]); // number of equations
  }


  
  // initialize parameters for problem
  PARAMETER param;
  param.ti = 0; // initial time
  param.tf = 10.0; // final time
  param.dt = (double)(param.tf - param.ti)/nt; // compute dt
  param.neq = neq; // number of equations
  param.nt = nt; // store number of time steps


  sol = new double[param.neq];
  // specify initial condition

  double xi;
  int Nx=param.neq/2;
  double dx = 1.0/(Nx+1);
  
  for (int i =0; i<Nx; i++) {
    xi = (i+1)*dx;
    sol[i]=1.0 + sin(2*3.14159265359*xi);
    sol[Nx+i] = 3.0;
  }


  // call ridc 
  ridc_be(order, param, sol);

  // output solution to screen
  for (int i = 0; i < param.neq; i++)
    printf("%14.12f\n", sol[i]);
  delete [] sol;
}
