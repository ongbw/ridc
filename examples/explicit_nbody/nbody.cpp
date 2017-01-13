#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <stdio.h>

#include "math.h"
#include "ridc.h"

using namespace std;

class ExplicitOde : public ODE {
public:
  ExplicitOde(int my_neq,
	      int my_nt,
	      double my_ti,
	      double my_tf,
	      double my_dt) {
    neq = my_neq;
    nt = my_nt;
    ti = my_ti;
    tf = my_tf;
    dt = my_dt;
  }
  
  void rhs(double t, double *u, double *f) {
    int npart = neq/2;
    int nion = npart/2;
    int nelec = npart/2;
    double d = 0.1;
    double q1 = 1.0/npart; // normalized ion charge
    double q2 = -1.0/npart; // normalized electron charge
    double m1 = 1000; // ions 1000x times heavier than electrons
    double m2 = 1;
    // assume the first npart/2 particles have mass m1 and charge q1

    // dx/dt = v
    for (int part =0; part<npart; part++) {
      f[part] = u[part + npart]; // velocity
    }

    int ind;
    double dist;
    // compute forces on ions
    for (int part = 0; part<nion; part++) {
      ind = part + npart;
      f[ind] = 0.0;
      for (int i = 0; i<nion; i++) {
	dist = u[part] - u[i];
	f[ind] = f[ind] + q1 / m1 * q1 * dist / sqrt(dist*dist + d*d);
      }
      for (int i = 0; i<nelec; i++) {
	dist = u[part] - u[nion+i];
	f[ind] = f[ind] + q1 / m1 * q2 * dist / sqrt(dist*dist + d*d);
      }    
    }
    //compute forces on electrons
    for (int part = 0; part<nelec; part++) {
      ind = part + nion + npart;
      f[ind] = 0.0;
      for (int i = 0; i<nion; i++) {
	dist = u[part+nion] - u[i];
	f[ind] = f[ind] + q2 / m2 * q1 * dist / sqrt(dist*dist + d*d);
      }
      for (int i = 0; i<nelec; i++) {
	dist = u[part+nion] - u[nion+i];
	f[ind] = f[ind] + q2 / m2 * q2 * dist / sqrt(dist*dist + d*d);
      }    
    }    
  }

  void step(double t, double * u, double * unew) {
    double* fold = new double[neq];
    rhs(t,u,fold);
  
    for (int i = 0; i < neq; i++)  {
      unew[i] = u[i] + dt*(fold[i]);
    }
    delete [] fold;
  }
};


int main(int argc, char *argv[]) {
  int order, nt;
  double *sol;

  if (argc != 3) {
    printf("usage: <executable> <order> <nt>  >  output_file\n");
    fflush(stdout);
    exit(1);
  } else {
    order = atoi(argv[1]); // order of method
    nt = atoi(argv[2]); // number of time steps
  }

  int npart = 4000;
  int neq = 2*npart;
  int ti = 0;
  int tf = 1;
  double dt = (double)(tf - ti)/nt; // compute dt
  
  // initialize ODE variable
  ExplicitOde *ode = new ExplicitOde(neq,nt,ti,tf,dt);

  int nion = npart/2;
  int nelec = npart/2;

  
  // specify initial condition
  sol = new double[neq];

  // initial ion location
  for (int i =0; i<nion; i++) {
    sol[i]=1.0/nion*(0.5+i);
  }
  // initial electron location
  for (int i =0; i<nelec; i++) {
    sol[nion+i]=sol[i];
  }
  // initial ion velocity
  for (int i =0; i<nion; i++) {
    sol[npart+i]=0.0;
  }
  // initial electron velocity
  for (int i =0; i<nelec; i++) {
    sol[npart+nion+i]=0.2*sin(2*M_PI*sol[nion+i]);
  }
  
  

  // call ridc 
  ridc_fe(ode,order, sol);

  // output solution to screen
  for (int i = 0; i < neq; i++)
    printf("%14.12f\n", sol[i]);
  delete [] sol;
  delete ode;
}
