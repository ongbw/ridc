#include "ode.h"
#include <stdio.h>

void init_params(const int& nt,
                 PARAMETER &param) {
    param.ti = 0;
    param.tf = 1;
    param.dt = (double)(param.tf - param.ti)/nt;
    param.neq = 2; // some neq
    param.nt = nt;
}

void initial_condition(PARAMETER param, double *u) {
  for (int i =0; i<param.neq; i++) {
    u[i]=1.0;
  }
}


void rhs(double t, double *u, PARAMETER param, double *f) {
  for (int i =0; i<param.neq; i++) {
      f[i]=-(i+1)*t*u[i];
  }
}

void step(double t, double* u, PARAMETER param, double* unew) {

  double* fold = new double[param.neq];
  rhs(t,u,param,fold);

  for (int i = 0; i < param.neq; i++)  {
    unew[i] = u[i] + param.dt*(fold[i]);
  }

  delete [] fold;
}
