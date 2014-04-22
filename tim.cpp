/***********************************************************************
*                   FAU Erlangen SS14
*                   Siwir2, Uebung 1 - Elliptic PDE with Multigrid
*                   written by Tim Brendel
*                   Januar 2014
************************************************************************/
#include "header.h"

using namespace std;


void save_in_file(const char *str, double *matrix, const int n_x, const int n_y, const double h_x, const double h_y){
	
	ofstream file;
	file.open(str, ios::out);
	if(!(file.is_open())){
		printf("%p konnte nicht gespeichert werden\n", str);
		exit(1);
	}
	//Sets the decimal precision to be used to format floating-point values on output operations.
	//New value for the decimal precision:12
	file << setprecision(12);
	for(int yi = 0; yi < n_y + 1; ++yi){
		for(int xj = 0; xj < n_x + 1; ++xj){
			file << xj * h_x<< '\t';
			file << yi * h_y<< '\t';
			file << matrix[yi * (n_x+1) + xj] << '\n';
		}
		file << endl;
	}
	file.close();
}

/*prolongation von grob/coarse nach fein/fine*/
void prolongation(double *u_co, double *u_fi, const int n_x, const int n_y){

	for(int yi = 0; yi < n_y; ++yi){
			for(int xj = 0; xj < n_x ; ++xj){
				if(yi > 0){
					if(xj > 0)
						u_fi[(yi - 1) * (n_x + 1) + xj - 1]	 += 1/4 * 1 * u_co[yi * (n_x + 1) + xj];
					u_fi[(yi - 1) * (n_x + 1) + xj] 		 += 1/4 * 2 * u_co[yi * (n_x + 1) + xj];
					u_fi[(yi - 1) * (n_x + 1) + xj + 1]	 += 1/4 * 1 * u_co[yi * (n_x + 1) + xj];
				}
				if(xj > 0){
					u_fi[yi * (n_x + 1) + xj - 1]			 += 1/4 * 2 * u_co[yi * (n_x + 1) + xj];
				}
				u_fi[yi * (n_x + 1) + xj] 				 += 1/4 * 4 * u_co[yi * (n_x + 1) + xj];
				u_fi[yi * (n_x + 1) + xj + 1]			 += 1/4 * 2 * u_co[yi * (n_x + 1) + xj];
				if(xj > 0){
					u_fi[(yi + 1) * (n_x + 1) + xj - 1]	 += 1/4 * 1 * u_co[yi * (n_x + 1) + xj];
				}
				u_fi[(yi + 1) * (n_x + 1) + xj]		 += 1/4 * 2 * u_co[yi * (n_x + 1) + xj];
				u_fi[(yi + 1) * (n_x + 1) + xj + 1]	 += 1/4 * 1 * u_co[yi * (n_x + 1) + xj];
			}
		}
}


void initialize_u_with_boundary_conditions(double *u, const int n_x, const int n_y, const double h_x, const double h_y){

	//beginne zu zaehlen links unten	
	//#pragma omp parallel for num_threads(16) if(n_x > 5000)
	for (int xj = 0; xj < n_x + 1; ++xj) { 
//		u[0 * (n_x + 1) + xj] = sin(PI * xj * h_x) * sinh( PI * 0 * h_y);	// == 0
		u[n_y * (n_x + 1) + xj] = sin(PI * xj * h_x) * sinh( PI * n_y * h_y);
	}
//	for(int yi = 0; yi < n_y + 1; ++yi){
//		u[yi * (n_x + 1) + 0] = sin(PI * 0 * h_x) * sinh(PI * yi * h_y);	// == 0
//		u[yi * (n_x + 1) + n_x] = sin(PI * n_x * h_x) * sinh(PI * yi * h_y);	// == 0
//	}

}


void do_gauss_seidel(double *u, double *f, const int n_x, const int n_y, const int c, 
						const double horizontal, const double vertical, const double center){

	/*do a gauss seidel iteration for c times
	 */
	for(int it = 0; it < c; it++ ){
		
// 		/*gauss seidel "normal"
// 		*/
// 		for(int yi = 1; yi < n_y; ++yi){
// 			for(int xj = 1; xj < n_x; ++xj){
// 				u[yi * (n_x + 1) + xj] = 	( f[yi * (n_x + 1) + xj]		+ horizontal * u[yi * (n_x + 1) + xj +1]
// 											+ horizontal * u[yi * (n_x + 1) + xj -1]
// 											+ vertical * u[(yi + 1) * (n_x + 1) + xj]
// 											+ vertical * u[(yi - 1) * (n_x + 1) + xj]
// 									) * center;
// 			}
// 		}
		/*red-black gauss seidel
		 */
		/*------red------*/
		//#pragma omp parallel for num_threads(32) schedule(static) firstprivate(u) if(n_x > 400 && n_y > 400)
		for(int yi = 1; yi < n_y ; yi++){
			for(int xj = 1 + (yi % 2); xj < n_x; xj += 2){
				u[yi * (n_x + 1) + xj] = 	( f[yi * (n_x + 1) + xj]	+ horizontal * u[yi * (n_x + 1) + xj +1]
										+ horizontal * u[yi * (n_x + 1) + xj -1]
										+ vertical * u[(yi + 1) * (n_x + 1) + xj]
										+ vertical * u[(yi - 1) * (n_x + 1) + xj]
								) * center;
			}
		}
		/*-------black-------*/
		//#pragma omp parallel for num_threads(32) schedule(static) firstprivate(u) if(n_x > 400 && n_y > 400)
		for(int yi = 1; yi < n_y; yi++) {
			for(int xj = 2 - (yi % 2); xj < n_x; xj += 2) {
				u[yi * (n_x + 1) + xj] = 	( f[yi * (n_x + 1) + xj] 	+ horizontal * u[yi * (n_x + 1) + xj +1]
										+ horizontal * u[yi * (n_x + 1) + xj -1]
										+ vertical * u[(yi + 1) * (n_x + 1) + xj]
										+ vertical * u[(yi - 1) * (n_x + 1) + xj]
								) * center;
			}
		}
	}
}


