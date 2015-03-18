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

  /// backward Euler update (simplistic is ODE is linear
  for (int i = 0; i < param.neq; i++)  {
    unew[i] = (1.0)/(1+param.dt*(t+param.dt)*(i+1))*u[i];
  }

  delete [] fold;
}
