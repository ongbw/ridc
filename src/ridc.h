/** 
    @file ridc.h
    @author Ong:Benjamin
    @version Revision 0.1
    @brief header file containing explanation of functions for the RIDC integrator
    @date 2015-03-18
*/


#ifndef _RIDC_H_
#define _RIDC_H_

#include <omp.h>
#include <cmath>
#include <algorithm>


struct PARAMETER {
  int neq;  /**< number of equations */
  int nt;  /**< number of time steps */
  double ti; /**< initial time */
  double tf;  /**< final time */
  double dt;  /**< time step */
};



template <class PARAMETER>
 void rhs(double t, double *u, PARAMETER param, double* f);
/**< This user-defined function (instantiated as a template) returns
   the right hand side of the ODE.  
   @return (by reference, f)
   @param param: structure containing number of equations
   @param t: current time step
   @param u: solution at current time
   @param f: returns f(t,y) by reference
*/

template <class PARAMETER>
void step(double t, double* u, PARAMETER param, double* unew);





  
// function declarations
void ridc_fe(int order, PARAMETER param, double *sol);
void ridc_be(int order, PARAMETER param, double *sol);
void lagrange_coeff(double *x, int Nx, int i, double *L);
double get_quad_weight(double *L, int Nx, double a, double b);
void integration_matrices(int N, double **S);
void init_unif_nodes(double *x, int Nx, double a, double b);
void corr_fe(double * uold,
             double ** fprev,
             double ** S,
             int index, int level,
             double t,
             PARAMETER param,
             double * unew);
void corr_be(double * uold,
             double ** fprev,
             double ** S,
             int index, int level,
             double t,
             PARAMETER param,
             double * unew);
#endif // _RIDC_H_
