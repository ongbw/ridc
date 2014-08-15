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
    int nx;
    int ny;
    double dt;
    double dx;
    double Lx;
    double Ly;
    double dy;
    double ti;
    double tf;
    double f1;
    double f2;
    double f3;
    double f4;
    double f5;
};

void init_params(const int& nt, PARAMETER &param);
void initial_condition(PARAMETER param, double* u);
double* rhs(double t, double *u, PARAMETER param);
void rhs(double t, double *u, PARAMETER param, double* f);
void newt(double t, double *uprev, double *uguess,
	  double *g, PARAMETER param);
void jac(double t, double *u, double *J, PARAMETER param);
void step(double t, double* u, PARAMETER param, double* unew);

#endif // _ODE_H_

