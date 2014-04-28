#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <math.h>
//#include <omp.h>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <cstring>


#define _USE_MATH_DEFINES

//Gauss-Seidel-Stencil
const double GS_HORIZONTAL=1.;
const double GS_VERTICAL=1.;
const double GS_CENTER=-4.;
//restriction Stencil
const double RES_HORIZONTAL=0.125;
const double RES_VERTICAL=0.125;
const double RES_CENTER=0.25;
const double RES_CORNER=0.0625;


#define LEADING_DIM n_x
#define IDX(i,j) ((j)*(LEADING_DIM)+(i))

static int l; // number of levels
static int n; // number of V-cycles
static double H; // meshsize
static int NX; // grid points in x-direction
static int NY; // grid points in y-direction


int getGridPointsDirichlet();
int getGridPointsNeumann();
void initializeGrid(double* u);
void initialize_u_with_boundary_conditions(double *u, const int n_x, const int n_y, const double h_x, const double h_y);
void initCoarseBD(const double* u_fi, double* u_co, int Nx_co);
void save_in_file(const char *str, double *matrix, const int n_x, const int n_y);

void do_gauss_seidel(double *u, double *f, const int n_x, const int n_y, const int c);
void residuum(double* res,double* f, double* u, const int n_x,const int n_y);
void restriction(double* f_co,double* res,const int n_x,const int n_y);
void mgm(double* u,double* f,int v1,int v2,int n_x, int n_y);
void prolongation(double *u_co, double *u_fi, const int n_x, const int n_y);
void calcResiduum(double *res, double *f, double *u, int n_x, int n_y);
void calculate_L2Norm(double *res, const double *u, const double *f, const int n_x, const int n_y);
double calcL2Norm(double *res, int n_x, int n_y);


