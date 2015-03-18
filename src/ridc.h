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
/**< This user-defined function (instantiated as a template) returns
   the Euler advance of the solution from time to to time t+dt
   @return (by reference, unew)
   @param param: structure containing number of equations, dt
   @param t: current time step
   @param u: solution at current time
   @param unew: returns unew by reference
*/


void ridc_fe(int order, PARAMETER param, double *sol);
/**< Main explicit ridc loop that initializes variables, integrates
   solution from ti to tf by bootstrapping the step function.
   @return (by reference) sol, the solution at the final time, param.tf
   @param param: structure containing number of equations, number of
   time steps, initial and final time, time step
   @param order: order of the RIDC method (predictor + number of correctors)
   @param sol: initial condition of the IVP
*/

void ridc_be(int order, PARAMETER param, double *sol);
/**< Main implicit ridc loop that initializes variables, integrates
   solution from ti to tf by bootstrapping the step function.

   @return (by reference) sol, the solution at the final time, param.tf

   @param param: structure containing number of equations, number of
   time steps, initial and final time, time step
   @param order: order of the RIDC method (predictor + number of correctors)
   @param sol: initial condition of the IVP
*/


void lagrange_coeff(double *x, int Nx, int i, double *L);
/**< RIDC helper function -- generates the coefficients for the
   lagrange interpolatory polynomials.

   @return (by reference) L: coefficients for the Lagrange
   intepolatory polynomial.  L is a vector of elements such that p(x)
   = L(0) + L(1)*x + L(2)*x^2 + ...
   @param x: quadrature nodes
   @param i: the i^{th} Lagrange polynomial
   @param Nx: number of quadrature nodes
   @param L: coefficients, returned by reference
*/


double get_quad_weight(double *L, int Nx, double a, double b);
/**< RIDC helper function -- generates quadrature weight,
   int(L_{n,i}(x),x=a..b)

   @return quadrature weights

   @param a: range of integration
   @param b: range of integration
   @param Nx: number of quadrature nodes
   @param L: coefficients for Lagrange poly,
	     L[0] + L[1]x + L[2]x^2 + ...
*/


void integration_matrices(int Nx, double **S);
/**< RIDC helper function -- constructions the integration matrix
   using get_quad_weight

   @return (by reference) the integration matrix S

   @param Nx: number of quadrature nodes
   @param S: integration matrix (by reference)

*/

void init_unif_nodes(double *x, int Nx, double a, double b);
/**< RIDC helper function -- initializes uniformly spaced quadrature nodes

   @return (by reference) x: uniformly spaced quadrature nodes

   @param Nx: number of quadrature nodes
   @param a: range of integration
   @param b: range of integration
   @param x: quadrature node location (returned by reference)

*/

void corr_fe(double * uold,
             double ** fprev,
             double ** S,
             int index, int level,
             double t,
             PARAMETER param,
             double * unew);
/**< RIDC helper function - solves error equation, updating the
   solution from time t to time t+param.dt.

   @return (by reference) unew: solution at time level t + param.dt

   @param uold: solution at time level t
   @param fprev: matrix containing derivative information from previous steps, previous level
   @param S: integration matrix (quadrature weights)
   @param index: decides which quadrature weights to use
   @param level: determines size of quadrature stencil
   @param t: current time iterate
   @param param: contains ODE parameters, for example, number of equations
   @param unew: solution at the new time level, passed by reference

*/

void corr_be(double * uold,
             double ** fprev,
             double ** S,
             int index, int level,
             double t,
             PARAMETER param,
             double * unew);
/**< RIDC helper function - solves error equation, updating the
   solution from time t to time t+param.dt.

   @return (by reference) unew: solution at time level t + param.dt

   @param uold: solution at time level t
   @param fprev: matrix containing derivative information from previous steps, previous level
   @param S: integration matrix (quadrature weights)
   @param index: decides which quadrature weights to use
   @param level: determines size of quadrature stencil
   @param t: current time iterate
   @param param: contains ODE parameters, for example, number of equations
   @param unew: solution at the new time level, passed by reference

*/

#endif // _RIDC_H_
