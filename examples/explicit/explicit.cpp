#include "ridc.h"
#include <stdio.h>

template <class PARAMTER>
void rhs(double t, double *u, PARAMETER param, double *f) {
  for (int i =0; i<param.neq; i++) {
      f[i]=-(i+1)*t*u[i];
  }
}

template <class PARAMETER>
void step(double t, double* u, PARAMETER param, double* unew) {

  double* fold = new double[param.neq];
  rhs(t,u,param,fold);

  for (int i = 0; i < param.neq; i++)  {
    unew[i] = u[i] + param.dt*(fold[i]);
  }

  delete [] fold;
}
