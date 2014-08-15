#include "ode.h"
#include <stdlib.h>
#include <stdio.h>
#include "mkl.h" // MKL Standard Lib
#include "mkl_lapacke.h" // MKL LAPACK
#include <cmath>

using namespace std;

void init_params(const int& nt,
                 PARAMETER &param) {
    param.ti = 0;
    param.tf = 10.0;
    param.dt = (double)(param.tf - param.ti)/nt;
    param.neq = 1e2; // some neq
    param.nt = nt;
}

void initial_condition(PARAMETER param, double *u) {
  double xi;
  int Nx=param.neq/2;
  double dx = 1.0/(Nx+1);

  for (int i =0; i<Nx; i++) {
    xi = (i+1)*dx;
    u[i]=1.0 + sin(2*3.14159265359*xi);
    u[Nx+i] = 3.0;
  }

}

void rhs(double t, double *u, PARAMETER param, double *f) {
  // brusselator

  double A = 1.0;
  double B = 3.0;
  double alpha = 0.02;
  double dx = 1.0/(param.neq+1);
  int Nx = param.neq/2;

  f[0] = A+u[0]*u[0]*u[Nx] -(B+1.0)*u[0] +alpha/dx/dx*(u[1]-2*u[0]+1.0);
  f[Nx-1] = A+u[Nx-1]*u[Nx-1]*u[2*Nx-1] -(B+1.0)*u[Nx-1]
    + alpha/dx/dx*(1.0-2*u[Nx-1]+u[Nx-2]);

  f[Nx] = B*u[0]- u[0]*u[0]*u[Nx] + alpha/dx/dx*(u[Nx+1]-2*u[Nx]+3.0);
  f[2*Nx-1] = B*u[Nx-1]- u[Nx-1]*u[Nx-1]*u[2*Nx-1]
    + alpha/dx/dx*(3.0-2*u[2*Nx-1]+u[2*Nx-2]);

  for (int i=1; i<Nx-1; i++) {
    f[i] = A + u[i]*u[i]*u[Nx+i] -(B+1.0)*u[i] +alpha/dx/dx*(u[i+1]-2*u[i]+u[i-1]);
    f[Nx+i] = B*u[i]- u[i]*u[i]*u[Nx+i] + alpha/dx/dx*(u[Nx+i+1]-2*u[Nx+i]+u[Nx+i-1]);
  }

}


///////////////////////////////////////////////////////
// Function Name: newt
// Usage: function for newton solve
//
///////////////////////////////////////////////////////
// Assumptions:
//
///////////////////////////////////////////////////////
// Inputs:
//         t, time
//         dt, timestep
//         uprev, solution at previous time level
//         uguess, solution for new time level
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//          g: function value at time t
//
///////////////////////////////////////////////////////
// Functions Called:
//
///////////////////////////////////////////////////////
// Called By:
//
///////////////////////////////////////////////////////
void newt(double t, double *uprev, double *uguess,
    double *g, PARAMETER param) {
  rhs(t,uguess,param,g);
  for (int i =0; i<param.neq; i++) {
    g[i] = uguess[i]-uprev[i]-param.dt*g[i];
  }
}


///////////////////////////////////////////////////////
// Function Name: jac
// Usage: compute approximation to jacobian
//
///////////////////////////////////////////////////////
// Assumptions:
//
///////////////////////////////////////////////////////
// Inputs:
//         t, independent variable
//         dt, timestep
//         u, dependent variable
//         param, problem parameters (d.g. neq, dt)
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//          J: approximate Jacobian at time t
//
///////////////////////////////////////////////////////
// Functions Called:
//
///////////////////////////////////////////////////////
// Called By: ???
//
///////////////////////////////////////////////////////
void jac(double t, double *u, double *J, PARAMETER param) {


  double d = 1e-5; // numerical jacobian approximation
  double *u1;
  double *f1;
  double *f;

  // need to allocate memory
  u1 = new double[param.neq];
  f1 = new double[param.neq];
  f = new double[param.neq];

  // note, instead of using newt, use rhs for efficiency
  rhs(t,u,param,f);

  for (int i = 0; i<param.neq; i++) {
    for (int j = 0; j<param.neq; j++) {
      u1[j] = u[j];
    }
    u1[i] = u1[i] + d;

    rhs(t,u1,param,f1);

    for (int j = 0; j<param.neq; j++) {
      J[j*param.neq+i] = -param.dt*(f1[j]-f[j])/d;
    }
    J[i*param.neq+i] = 1.0+J[i*param.neq+i];
  }

  // need to delete memory

  delete [] u1;
  delete [] f1;
  delete [] f;

}


void step(double t, double* uold, PARAMETER param, double* unew) {

  // t is current time
  double tnew = t + param.dt;

  double NEWTON_TOL = 1.0e-14;
  int NEWTON_MAXSTEP = 1000;

  double* stepsize = new double[param.neq];
  double* uguess = new double[param.neq];

  double *J = new double[param.neq*param.neq];

  for (int i=0; i<param.neq; i++) {
    uguess[i] = uold[i];
  }

  double maxstepsize;

  int counter=0;


  while (1) {

    // create rhs
    newt(tnew,uold,uguess,stepsize,param);

    // compute numerical Jacobian
    jac(tnew,uguess,J,param);

    int info;
    int n = param.neq;
    int nrhs = 1;

    int pivot[param.neq];

    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, param.neq, nrhs, J, n, pivot, stepsize, 1);

    // check for convergence
    maxstepsize = 0.0;
    for (int i = 0; i<param.neq; i++) {
      uguess[i] = uguess[i] - stepsize[i];
      if (abs(stepsize[i])>maxstepsize) {
	maxstepsize = abs(stepsize[i]);
      }
    }

    // if update sufficiently small enough
    if ( maxstepsize < NEWTON_TOL) {
      break;
    }

    counter++;
    //error if too many steps
    if (counter > NEWTON_MAXSTEP) {
      fprintf(stderr,"max newton iterations reached\n");
      exit(42);
    }
  } // end newton iteration

  for (int i=0; i<param.neq; i++) {
    unew[i] = uguess[i];
  }


  delete [] uguess;

  delete [] stepsize;

  delete [] J;
}
