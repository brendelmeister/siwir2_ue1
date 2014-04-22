/***********************************************************************
*                   FAU Erlangen SS14
*                   Siwir2, Uebung 1 - Elliptic PDE with Multigrid
*                   written by Tim Brendel
*                   Januar 2014
************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <sys/time.h>
#include <iostream>


#define PI M_PI
/*TODO: globale Variablen: NX, NY, gamma, l, n, h*/

void save_in_file(const char *str, double *matrix, const int n_x, const int n_y, const double h_x, const double h_y);
void prolongation(double *u_co, double *u_fi, const int n_x, const int n_y);
void initialize_u_with_boundary_conditions(double *u, const int n_x, const int n_y, const double h_x, const double h_y);
void do_gauss_seidel(double *u, double *f, const int n_x, const int n_y, const int c, 
						const double horizontal, const double vertical, const double center);
int do_nothing(void);
void do_something_extremely_useful(void);
