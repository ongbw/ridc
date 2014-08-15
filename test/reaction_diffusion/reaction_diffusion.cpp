#include "ode.h"
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

using namespace std;

void init_params(const int& nt,
                 PARAMETER &param) {
    param.ti = 0;
    param.tf = 0.75;
    param.dt = (double)(param.tf - param.ti)/nt;
    param.nx = 50;
    param.Lx = 1.0;
    param.ny = 50;
    param.Ly = 1.0;
    param.dx = (double) param.Lx/(param.nx-1);
    param.dy = (double) param.Ly/(param.ny-1);
    param.neq = 2*param.nx*param.ny; // 2 times the total number of nodes
    param.nt = nt;
// parameters for f
    param.f1 = 0.1305; //a 
    param.f2 = 0.7695; //%b
    param.f3 = 100;  //%k
    param.f4=0.05; //%D1
    param.f5 =1;//%D2

} // end of init_params



void initial_condition(PARAMETER param, double *u) {
  int nyy = param.ny;
  int nxx = param.nx;
   double dx = param.dx;
   double dy = param.dy;
 

// new stuff

	for (int j=0; j < nyy; j++) {
		for (int i=0; i < nxx; i++) {
        		u[j*nxx+i] = param.f1+param.f2+(1e-3)*exp(-100*(pow((i*dx-1/3),2)+pow((j*dy-1/2),2)));
        		u[j*nxx+i+nxx*nyy] = param.f2/pow((param.f1+param.f2),2);
		}
	}
} // end of initial condition function 




void rhs(double t, double *w, PARAMETER param, double *f) {
	// 2D reaction diffusion system

	int Nx = param.nx;
	int Ny = param.ny;
	double a = param.f1; 
        double b= param.f2; 
        double k=param.f3;
         double D1=param.f4; 
        double D2=param.f5;
        double dx=param.dx; 
        double dy=param.dy; 
	double dx2 = pow(dx,2.0);
	double dy2 = pow(dy,2.0);
	int i;
	int j;


//% take care of interior nodes

for (int j = 1; j < Ny-1; j++)
    {
    for (int i =1; i < Nx-1; i++)
        {
        f[(j)*Nx+i] = D1/dx2*( w[(j)*Nx+i-1]-2*w[(j)*Nx+i]+w[(j)*Nx+i+1]) + D1/dy2*(w[(j-1)*Nx+i]-2*w[(j)*Nx+i]+w[(j+1)*Nx+i]) +k*(a-w[(j)*Nx+i]+pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);
        f[(j)*Nx+i+Nx*Ny] = D2/dx2*(w[(j)*Nx+i-1+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j)*Nx+i+1+Nx*Ny]) + D2/dy2*(w[(j-1)*Nx+i+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j+1)*Nx+i+Nx*Ny]) +k*(b-pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);
        }
    }

// interior done

    // take care of the edges
    
    // bottom edge
	j=0;
    for (int i=1; i< Nx-1;i++)
        {
 f[(j)*Nx+i] = D1/dx2*( w[(j)*Nx+i-1]-2*w[(j)*Nx+i]+w[(j)*Nx+i+1]) + D1/dy2*(w[1*Nx+i]-2*w[(j)*Nx+i]+w[(j+1)*Nx+i]) +k*(a-w[(j)*Nx+i]+pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);
        f[(j)*Nx+i+Nx*Ny] = D2/dx2*(w[(j)*Nx+i-1+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j)*Nx+i+1+Nx*Ny]) + D2/dy2*(w[1*Nx+i+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j+1)*Nx+i+Nx*Ny]) +k*(b-pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);
        
        }
   // done 
    
   // bottom left corner 
	j=0;
	i=0;

 f[(j)*Nx+i] = D1/dx2*( w[(j)*Nx+1]-2*w[(j)*Nx+i]+w[(j)*Nx+i+1]) + D1/dy2*(w[1*Nx+i]-2*w[(j)*Nx+i]+w[(j+1)*Nx+i]) +k*(a-w[(j)*Nx+i]+pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);
        f[(j)*Nx+i+Nx*Ny] = D2/dx2*(w[(j)*Nx+1+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j)*Nx+i+1+Nx*Ny]) + D2/dy2*(w[1*Nx+i+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j+1)*Nx+i+Nx*Ny]) +k*(b-pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);

   // done 

 // bottom right corner
	j=0;
	i=Nx-1;
 f[(j)*Nx+i] = D1/dx2*( w[(j)*Nx+i-1]-2*w[(j)*Nx+i]+w[(j)*Nx+Nx-2]) + D1/dy2*(w[1*Nx+i]-2*w[(j)*Nx+i]+w[(j+1)*Nx+i]) +k*(a-w[(j)*Nx+i]+pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);
        f[(j)*Nx+i+Nx*Ny] = D2/dx2*(w[(j)*Nx+i-1+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j)*Nx+Nx-2+Nx*Ny]) + D2/dy2*(w[1*Nx+i+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j+1)*Nx+i+Nx*Ny]) +k*(b-pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);


  // done 
    
    // top edge
    
	j=Ny-1;
    
    for (int i=1; i < Nx-1; i++)
        {

 f[(j)*Nx+i] = D1/dx2*( w[(j)*Nx+i-1]-2*w[(j)*Nx+i]+w[(j)*Nx+i+1]) + D1/dy2*(w[(j-1)*Nx+i]-2*w[(j)*Nx+i]+w[(Ny-2)*Nx+i]) +k*(a-w[(j)*Nx+i]+pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);
        f[(j)*Nx+i+Nx*Ny] = D2/dx2*(w[(j)*Nx+i-1+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j)*Nx+i+1+Nx*Ny]) + D2/dy2*(w[(j-1)*Nx+i+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(Ny-2)*Nx+i+Nx*Ny]) +k*(b-pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);

        }

   // done  
// top left corner
	j=Ny-1;
	i=0;

 f[(j)*Nx+i] = D1/dx2*( w[(j)*Nx+1]-2*w[(j)*Nx+i]+w[(j)*Nx+i+1]) + D1/dy2*(w[(j-1)*Nx+i]-2*w[(j)*Nx+i]+w[(Ny-2)*Nx+i]) +k*(a-w[(j)*Nx+i]+pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);
        f[(j)*Nx+i+Nx*Ny] = D2/dx2*(w[(j)*Nx+1+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j)*Nx+i+1+Nx*Ny]) + D2/dy2*(w[(j-1)*Nx+i+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(Ny-2)*Nx+i+Nx*Ny]) +k*(b-pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);

    
// done
   // top right corner 
	j=Ny-1;
	i=Nx-1;
    
 f[(j)*Nx+i] = D1/dx2*( w[(j)*Nx+i-1]-2*w[(j)*Nx+i]+w[(j)*Nx+Nx-2]) + D1/dy2*(w[(j-1)*Nx+i]-2*w[(j)*Nx+i]+w[(Ny-2)*Nx+i]) +k*(a-w[(j)*Nx+i]+pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);
        f[(j)*Nx+i+Nx*Ny] = D2/dx2*(w[(j)*Nx+i-1+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j)*Nx+Nx-2+Nx*Ny]) + D2/dy2*(w[(j-1)*Nx+i+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(Ny-2)*Nx+i+Nx*Ny]) +k*(b-pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);

    
   // done 
    
    // left edge
    
	i=0;
    
    for (int j=1; j < Ny-1; j++)
        {
 f[(j)*Nx+i] = D1/dx2*( w[(j)*Nx+1]-2*w[(j)*Nx+i]+w[(j)*Nx+i+1]) + D1/dy2*(w[(j-1)*Nx+i]-2*w[(j)*Nx+i]+w[(j+1)*Nx+i]) +k*(a-w[(j)*Nx+i]+pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);
        f[(j)*Nx+i+Nx*Ny] = D2/dx2*(w[(j)*Nx+1+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j)*Nx+i+1+Nx*Ny]) + D2/dy2*(w[(j-1)*Nx+i+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j+1)*Nx+i+Nx*Ny]) +k*(b-pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);


        }
// done
    
	// bottom left corner again?
	i=0;
	j=0;
    f[(j)*Nx+i] = D1/dx2*( w[(j)*Nx+1]-2*w[(j)*Nx+i]+w[(j)*Nx+i+1]) + D1/dy2*(w[(1)*Nx+i]-2*w[(j)*Nx+i]+w[(j+1)*Nx+i]) +k*(a-w[(j)*Nx+i]+pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);
        f[(j)*Nx+i+Nx*Ny] = D2/dx2*(w[(j)*Nx+1+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j)*Nx+i+1+Nx*Ny]) + D2/dy2*(w[(1)*Nx+i+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j+1)*Nx+i+Nx*Ny]) +k*(b-pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);

 
   // done 
    
    
// top left corner again?
	i=0;
	j=Ny-1;
    
 f[(j)*Nx+i] = D1/dx2*( w[(j)*Nx+1]-2*w[(j)*Nx+i]+w[(j)*Nx+i+1]) + D1/dy2*(w[(j-1)*Nx+i]-2*w[(j)*Nx+i]+w[(Ny-2)*Nx+i]) +k*(a-w[(j)*Nx+i]+pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);
        f[(j)*Nx+i+Nx*Ny] = D2/dx2*(w[(j)*Nx+1+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j)*Nx+i+1+Nx*Ny]) + D2/dy2*(w[(j-1)*Nx+i+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(Ny-2)*Nx+i+Nx*Ny]) +k*(b-pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);

// done
    
    
    // right edge
	i=Nx-1;
    
     for (int j=1; j< Ny-1;j++)
         {

 f[(j)*Nx+i] = D1/dx2*( w[(j)*Nx+i-1]-2*w[(j)*Nx+i]+w[(j)*Nx+Nx-2]) + D1/dy2*(w[(j-1)*Nx+i]-2*w[(j)*Nx+i]+w[(j+1)*Nx+i]) +k*(a-w[(j)*Nx+i]+pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);
        f[(j)*Nx+i+Nx*Ny] = D2/dx2*(w[(j)*Nx+i-1+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j)*Nx+Nx-2+Nx*Ny]) + D2/dy2*(w[(j-1)*Nx+i+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j+1)*Nx+i+Nx*Ny]) +k*(b-pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);


     }   
// done
    
// right bottom corner
	i=Nx-1;
	j=0;
    
 f[(j)*Nx+i] = D1/dx2*( w[(j)*Nx+i-1]-2*w[(j)*Nx+i]+w[(j)*Nx+Nx-2]) + D1/dy2*(w[(1)*Nx+i]-2*w[(j)*Nx+i]+w[(j+1)*Nx+i]) +k*(a-w[(j)*Nx+i]+pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);
        f[(j)*Nx+i+Nx*Ny] = D2/dx2*(w[(j)*Nx+i-1+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j)*Nx+Nx-2+Nx*Ny]) + D2/dy2*(w[(1)*Nx+i+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j+1)*Nx+i+Nx*Ny]) +k*(b-pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);

// done
    
   // right top corner
	i=Nx-1;
	j=Ny-1;
    
 f[(j)*Nx+i] = D1/dx2*( w[(j)*Nx+i-1]-2*w[(j)*Nx+i]+w[(j)*Nx+Nx-2]) + D1/dy2*(w[(j-1)*Nx+i]-2*w[(j)*Nx+i]+w[(Ny-2)*Nx+i]) +k*(a-w[(j)*Nx+i]+pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);
        f[(j)*Nx+i+Nx*Ny] = D2/dx2*(w[(j)*Nx+i-1+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(j)*Nx+Nx-2+Nx*Ny]) + D2/dy2*(w[(j-1)*Nx+i+Nx*Ny]-2*w[(j)*Nx+i+Nx*Ny]+w[(Ny-2)*Nx+i+Nx*Ny]) +k*(b-pow(w[(j)*Nx+i],2)*w[(j)*Nx+i+Nx*Ny]);


   // done 
    
}    // end of f function for right hand side 
    
    


void step(double t, double* u, PARAMETER param, double* unew) {

  // initialize solution
	double fold[param.neq];
 rhs(t, u, param,fold);

  for (int i = 0; i < param.neq; i++)  {
    unew[i] = u[i] + param.dt*(fold[i]);
  }

//  delete [] fold;
}
