/**
 * 2012-11-09
 *
 * ridc.h
 */

#ifndef _RIDC_H_
#define _RIDC_H_

#include <omp.h>
#include <cmath>
#include <algorithm>

#include "ode.h"

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
