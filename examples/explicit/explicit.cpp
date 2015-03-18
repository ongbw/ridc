#include "ridc.h"
#include <stdio.h>

template <>
void rhs(double t, double *u, PARAMETER param, double *f) {
  /// return rhs of the odes d/dt(u_i)  = -(i+1) t u_i
  for (int i =0; i<param.neq; i++) {
      f[i]=-(i+1)*t*u[i];
  }
}

template <>
void step(double t, double* u, PARAMETER param, double* unew) {

  double* fold = new double[param.neq];
  rhs(t,u,param,fold);

  /// forward euler update
  for (int i = 0; i < param.neq; i++)  {
    unew[i] = u[i] + param.dt*(fold[i]);
  }

  delete [] fold;
}
