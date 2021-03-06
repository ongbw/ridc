#include "ridc.h"
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <gsl/gsl_linalg.h>

#ifndef _BRUSSELATOR_H_
#define _BRUSSELATOR_H_

class Brusselator_GSL : public ODE {
public:
  Brusselator_GSL(int my_neq, int my_nt, double my_ti, double my_tf, double my_dt) {
    neq = my_neq; 
    nt = my_nt;  
    ti = my_ti; 
    tf = my_tf; 
    dt = my_dt; 
  }
  
  void rhs(double t, double *u, double *f) {
    /** user implemented rhs function, u'=rhs(t,u) 
       @return (by reference) f: rhs(t,u)
       @param t current time step
       @param u solution u at time t
       @param f rhs(t,u)
    */    
    double A = 1.0;
    double B = 3.0;
    double alpha = 0.02;
    double dx = 1.0/(neq+1);
    int Nx = neq/2;
    
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

  void step(double t, double * u, double * unew) {
    /** user implemented step function, for advancing the solution from t to t+dt 
       @return (by reference) unew: solution at time t+dt
       @param t current time step
       @param u solution u at time t
       @param unew solution at time t+dt
    */
    double tnew = t + dt;

    double NEWTON_TOL = 1.0e-14;
    int NEWTON_MAXSTEP = 1000;

    double * J;
    double * stepsize;
    double * uguess;
    
    stepsize = new double[neq];
    uguess = new double[neq];

    J = new double[neq*neq];

    // GSL
    gsl_matrix_view m;
    gsl_vector_view b;
    gsl_vector *x = gsl_vector_alloc(neq);
    gsl_permutation *p = gsl_permutation_alloc(neq);

  
    // set initial condition
    for (int i=0; i<neq; i++) {
      uguess[i] = u[i];
    }

    double maxstepsize;
    
    int counter=0;
    int * pivot = new int[neq];

    while (1) {

      // create rhs
      newt(tnew,u,uguess,stepsize);

      // compute numerical Jacobian
      jac(tnew,uguess,J);

      m = gsl_matrix_view_array(J,neq,neq);
      b = gsl_vector_view_array(stepsize,neq);
      
      int s; // presumably, this is some flag info
      
      gsl_linalg_LU_decomp (&m.matrix,p,&s);
      gsl_linalg_LU_solve (&m.matrix,p,&b.vector, x);
      
      
      // check for convergence
      maxstepsize = 0.0;
      for (int i = 0; i<neq; i++) {
	uguess[i] = uguess[i] - x->data[i];
	if (std::abs(stepsize[i])>maxstepsize) {
	  maxstepsize = std::abs(x->data[i]);
	}
      }
      
      // if update sufficiently small enough
      if ( maxstepsize < NEWTON_TOL) {
	break;
      }

      counter++;
      //error if too many steps
      if (counter > NEWTON_MAXSTEP) {
	fprintf(stderr,"max Newton iterations reached\n");
	exit(42);
      }
    } // end Newton iteration
    
    for (int i=0; i<neq; i++) {
      unew[i] = uguess[i];
    }

    delete [] pivot;
    delete [] uguess;
    delete [] stepsize;
    delete [] J;
    
    gsl_permutation_free (p);
    gsl_vector_free (x);
  }

 private:
  void newt(double t, double *uprev, double *uguess,
	    double *g){
    /** Helper function to compute the next Newton step
       for solving a system of equations
       
       @return (by reference) g how far from zero we are
       @param t current time step
       @param uguess current solution guess
       @param uprev solution at previous time step
       @param g how far from zero we are, returned by reference
    */
    rhs(t,uguess,g);
    for (int i =0; i<neq; i++) {
      g[i] = uguess[i]-uprev[i]-dt*g[i];
    }
  }
    
  void jac(double t, double *u, double *J){
    /** Helper function to the Jacobian matrix (using finite differences)
       for advancing the solution from time t(n) to t(n+1) using an
       implicit Euler step on a system of equations
       
       @return (by reference) J the Jacobian for the Newton step
       @param t current time step
       @param u function value at the current time step
       @param J Jacobian, returned by reference
    */
    double d = 1e-5; // numerical Jacobian approximation
    double *u1;
    double *f1;
    double *f;
    
    
    // need to allocate memory
    u1 = new double[neq];
    f1 = new double[neq];
    f = new double[neq];
    
    // note, instead of using newt, use rhs for efficiency
    rhs(t,u,f);
    
    for (int i = 0; i<neq; i++) {
      for (int j = 0; j<neq; j++) {
	u1[j] = u[j];
      }
      u1[i] = u1[i] + d;
      
      rhs(t,u1,f1);
      
      for (int j = 0; j<neq; j++) {
	J[j*neq+i] = -dt*(f1[j]-f[j])/d;
      }
      J[i*neq+i] = 1.0+J[i*neq+i];
    }
    
    // need to delete memory
    delete [] u1;
    delete [] f1;
    delete [] f;
  }
  
};

#endif // _BRUSSELATOR_H_
