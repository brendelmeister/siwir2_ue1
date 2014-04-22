#include <stdlib.h>
#include <iostream>
#include <math.h>

using namespace std;
static int l; // number of levels
static int n; // number of V-cycles
static double h; // meshsize
static int NX; // grid points in x-direction
static int NY;  // grid points in y-direction
static int gam = 2; //cycle-structure

int getGridPointsDirichlet();
void initializeGrid(double* u);
void RBGS(double* u); // Red Black Gauss Seidel Relaxation
void bilinearInterpolation(double *u, int coarseLevel, int fineLevel); // Prolongation
void fullWeighting(double *u); // Restriction

int main(int argc, char *argv[]) {
	if (argc != 3) {
		cout<<"Wrong number of arguments!"<<endl;
		exit(EXIT_FAILURE);
	}
	l = atoi(argv[1]);
	n = atoi(argv[2]);
	NX = NY = getGridPointsDirichlet();
	h = 1/NX;

	double* u = new double[NX*NY]; // grid
	initializeGrid(u);
}

int getGridPointsDirichlet(){
	return pow(2,l)+1;
}

void initializeGrid(double* u){
	for(int i=0 ;i<NX;++i){
		for (int j = 0; j < NY; ++j) {
			u[i*NY+j] = 0; 
			if(i == NX-1 || j == NY-1) // as sin = 0 for i==0 or j == 0
			{
				u[i*NY+j] = sin(M_PI*i*h)*sinh(M_PI*j*h);
			}
		}
	}
}

void RBGS(double* u){

}

void bilinearInterpolation(double *u, int coarseLevel, int fineLevel){

}

void calcResiduum(double *res, double *u, double *f){

	for(int i= 0; i<NX; i++){
		for(int j =0 ; j<NY; j++){
			res[i*NY+j] = f[i*NY+j] - ((1.0/(h*h))*(-4.0*u[i*NY+j]+u[(i+1)*NY+j]+u[(i-1)*NY+j]+u[i*NY+j+1]+u[i*NY+j-1]));
		}
	}
}

void calcL2Norm(double *res){
	double norm = 0;
	for(int i= 0; i<NX; i++){
		for(int j =0 ; j<NY; j++){
			norm += res[i*NY+j]*res[i*NY+j];
		}
	}
	norm = sqrt(norm)/(NX*NY);	
}
/*void measureError(double *u, double *u_co, double gridsize){
  double *error = new double[NX*NY];
  for(int i = 0; i<NX; i++){
  for(int j = 0; j < NY; j++){
  error[i*NY+NX] = sqrt(u[i*NY+j]*u[i*NX+j]-u_co[i*NY+j]*u_co[i*NY+j]);
  }
  }
  }*/
