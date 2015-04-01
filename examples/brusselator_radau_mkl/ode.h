/**
 * 2014--6--18
 * ode.h
 */
#ifndef _ODE_H_
#define _ODE_H_

/**
 * Requires at least
 * dt - delta t
 * neq - number of equations
 */
struct PARAMETER {
    int neq;
    int nt;
    double ti;
    double tf;
    double dt;
};

struct BUTCHER {
  int S;
  double * b;
  double * c;
  double ** A;
  
};

void rhs(double t, double *u, PARAMETER param, double* f);
void newt(double t, double *uprev, double * Kguess,
	  double *g, PARAMETER param, BUTCHER rk);
void jac(double t, double *uprev, double * Kguess,
	 double *J, PARAMETER param, BUTCHER rk);
void step(double t, double* u, PARAMETER param, double* unew, BUTCHER rk);

#endif // _ODE_H_
