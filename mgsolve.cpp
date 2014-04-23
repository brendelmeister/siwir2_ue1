/****************************************************************************
*                   FAU Erlangen SS14
*                   Siwir2 Uebung 1 - Elliptic PDE with Multigrid
*                   written by Lina Gundelwein, Michael Hildebrand, Tim Brendel
*                   April 2014
*****************************************************************************/

#include "header.h"
using namespace std;


int main(int argc, char *argv[]) {
	if (argc != 3) {
		cout<<"Wrong number of arguments!"<<endl;
		cout<<"call ./mgsolve l n "<<endl;
		cout<<"l: number of levels; n: numer of V-cycles"<<endl;
		exit(EXIT_FAILURE);
	}
	l = atoi(argv[1]);
	n = atoi(argv[2]);
	NX = NY = getGridPointsDirichlet();
	H = 1.0/(NX-1);


	double* u = new double[NX*NY]; // initialise arrays
	initializeGrid(u);
	double* f = new double[NX*NY];
	memset(f,0,sizeof(double)*NY*NX);
	do_gauss_seidel(u,f,NX,NY,25);
	
	
// 	for(int i=0;i<n;i++){ //multigrid steps
// 		mgm( u, f,2,1,NX, NY);
// 	}
	save_in_file("boundaries.txt", u, NX, NY);
    delete[] u;
    delete[] f;
}

void save_in_file(const char *str, double *matrix, const int n_x, const int n_y){
	ofstream file;
	file.open(str, ios::out);
	if(!(file.is_open())){
		printf("%p konnte nicht gespeichert werden\n", str);
		exit(1);
	}
	double hx_local = 1.0/(n_x-1);
	double hy_local = 1.0/(n_y-1);
	//Sets the decimal precision to be used to format floating-point values on output operations.
	//New value for the decimal precision:12
	file << setprecision(12);
	for(int yi = 0; yi < n_y ; ++yi){
		for(int xj = 0; xj < n_x ; ++xj){
			file << xj*hx_local << '\t';
			file << yi*hy_local << '\t';
			file << matrix[yi * (n_x) + xj] << '\n';
		}
		file << endl;
	}
	file.close();
}

/*prolongation von grob/coarse nach fein/fine*/
void prolongation(double *u_co, double *u_fi, const int n_x, const int n_y){

	for(int yi = 0; yi < n_y; ++yi){
		for(int xj = 0; xj < n_x ; ++xj){
			if(yi>0){
				if(xj > 0){
					// south west
					u_fi[(yi - 1) * (n_x) + xj - 1]	= 1/4 * 1 * u_co[yi * (n_x) + xj];
				}	
				if(xj<n_x-1){	
					// south east
					u_fi[(yi - 1) * (n_x) + xj + 1]	= 1/4 * 1 * u_co[yi * (n_x) + xj];	
				}
				// south
				u_fi[(yi - 1) * (n_x) + xj] = 1/4 * 2 * u_co[yi * (n_x) + xj];

			}
			if(xj>0){
				// west
				u_fi[yi * (n_x) + xj - 1] = 1/4 * 2 * u_co[yi * (n_x) + xj];
				if(yi < n_y-1){	
					// north west
					u_fi[(yi + 1) * (n_x) + xj - 1]	= 1/4 * 1 * u_co[yi * (n_x) + xj];
				}
			}
			// centre
			u_fi[yi * (n_x) + xj]= 1/4 * 4 * u_co[yi * (n_x) + xj];

			if(xj<n_x-1){
				// east
				u_fi[yi * (n_x) + xj + 1]= 1/4 * 2 * u_co[yi * (n_x) + xj];
			}
			if(yi < n_y-1){
				// north
				u_fi[(yi + 1) * (n_x) + xj]= 1/4 * 2 * u_co[yi * (n_x) + xj];
				// north east
				u_fi[(yi + 1) * (n_x) + xj + 1]	= 1/4 * 1 * u_co[yi * (n_x) + xj];
			}
		}
	}
}


void do_gauss_seidel(double *u, double *f, const int n_x, const int n_y, const int c){

	if(n_x != n_y){
		cout<<"error: grid not quadratic"<<endl;
		exit(EXIT_FAILURE);
	}
	double h = 1.0 / n_x;
/*

	//do a gauss seidel iteration for c times
	for(int it = 0; it < c; ++it ){

	 // gauss seidel "normal" 
	 for(int yi = 1; yi < n_y; ++yi){
	  for(int xj = 1; xj < n_x; ++xj){
	   u[yi * n_x + xj] = ( f[yi * n_x + xj]
	     + GS_HORIZONTAL * u[yi * n_x + xj +1]
	     + GS_HORIZONTAL * u[yi * n_x + xj -1]
	     + GS_VERTICAL * u[(yi + 1) * n_x + xj]
	     + GS_VERTICAL * u[(yi - 1) * n_x + xj]
	     ) * GS_CENTER / (h*h);
	  }
	 }
	 // 		//red-black gauss seidel 
	 // 		//------red------
	 // 		//#pragma omp parallel for num_threads(32) schedule(static) firstprivate(u) if(n_x > 400 && n_y > 400)
	 // 		for(int yi = 1; yi < n_y ; yi++){
	 // 			for(int xj = 1 + (yi % 2); xj < n_x; xj += 2){
	 // 				u[yi * n_x + xj] = ( f[yi * n_x + xj]	+ GS_HORIZONTAL * u[yi * n_x + xj +1]
	 // 																		+ GS_HORIZONTAL * u[yi * n_x + xj -1]
	 // 																		+ GS_VERTICAL * u[(yi + 1) * n_x + xj]
	 // 																		+ GS_VERTICAL * u[(yi - 1) * n_x + xj]
	 // 											) * GS_CENTER / (h*h);
	 // 			}
	 // 		}
	 // 		//-------black-------
	 // 		//#pragma omp parallel for num_threads(32) schedule(static) firstprivate(u) if(n_x > 400 && n_y > 400)
	 // 		for(int yi = 1; yi < n_y; yi++) {
	 // 			for(int xj = 2 - (yi % 2); xj < n_x; xj += 2) {
	 // 				u[yi * (n_x + 1) + xj] = ( f[yi * n_x + xj] 	+ GS_HORIZONTAL * u[yi * n_x + xj +1]
	 // 																				+ GS_HORIZONTAL * u[yi * n_x + xj -1]
	 // 																				+ GS_VERTICAL * u[(yi + 1) * n_x + xj]
	 // 																				+ GS_VERTICAL * u[(yi - 1) * n_x + xj]
	 // 													) * GS_CENTER / (h*h);
	 // 			}
	 // 		}
	}
	*/

	//#pragma omp parallel
	for (int iter=0;iter<c;iter++)
	{
	 //red
//#pragma omp for schedule(static)
	 for (int y=1;y<n_y-1;y++)
	 {
	  for (int x=(y%2)+1;x<n_x-1;x+=2)
	  {
	   u[IDX(x,y)] = 1.0/4.0 * (h*h*f[IDX(x,y)] + (u[IDX(x,y-1)] + u[IDX(x,y+1)] + u[IDX(x-1,y)] + u[IDX(x+1,y)]));
	  }
	 }
	 //black
//#pragma omp for schedule(static)
	 for (int y=1;y<n_y-1;y++)
	 {
	  for (int x=((y+1)%2)+1;x<n_x-1;x+=2)
	  {
	   u[IDX(x,y)] = 1.0/4.0 * (h*h*f[IDX(x,y)] + (u[IDX(x,y-1)] + u[IDX(x,y+1)] + u[IDX(x-1,y)] + u[IDX(x+1,y)]));
	  }
	 }

	}
}


int getGridPointsDirichlet(){
	return pow(2,l)+1;
}

void initializeGrid(double* u){
// 	for(int j=0 ;j<NY;++j){
// 		for (int i = 0; i < NX; ++i) {
// 			if(i == NX-1) // as sin = 0 for i==0 or j == 0
// 			{
// 				u[j*NX+i] = sin(M_PI*j*h)*sinh(M_PI*i*h);
// 			}
// 		}
// 	}
	for(int i = 0; i < NX; ++i){
		u[(NY-1) * NX + i] = sin(M_PI*i*H) * sinh(M_PI*(NY-1)*H);
	}
}


//do restriction from residuum to f_coarse
void restriction(double* f_co,double* res,const int n_x,const int n_y,
		const double horizontal, const double vertical, const double center, const double corner){
	int Nx_co=(n_x/2)+1;
	int Ny_co=(n_y/2)+1;
	for(int j=0;j<Ny_co;j++){
		for(int i=0;i<Nx_co;i++){
			f_co[j*Ny_co+i] =RES_CENTER*res[j*2*n_x+i*2]+ //restriction stencil
				RES_HORIZONTAL*(res[(j*2*n_x+i*2)-1]+ res[(j*2*n_x+i*2)+1])+
				RES_VERTICAL*(res[((j*2-1)*n_x+i*2)]+ res[((j*2+1)*n_x+i*2)])+
				RES_CORNER*(res[((j*2-1)*n_x+i*2)-1]+ res[((j*2-1)*n_x+i*2)+1]+
						res[((j*2+1)*n_x+i*2)-1]+ res[((j*2+1)*n_x+i*2)+1]);
		}
	}
}

// copy boundaries from fine to coarse grid
void initCoarseBD(const double* u_fi, double* u_co, int Nx_co){
	for(int j=0;j<Nx_co;j++){
		u_co[j]= u_fi[2*j];
	}

}


//recursive multigrid function
void mgm(double* u,double* f,int v1,int v2,int n_x, int n_y){

	do_gauss_seidel(u,f,n_x,n_y,v1);//Pre-smoothing

	double* res = new double[n_y*n_x];

	residuum(res, f, u, n_x, n_y); //residuum calculation

	int Nx_co=(n_x/2)+1; //calculating coarse grid size
	int Ny_co=(n_y/2)+1;
	double* f_co=new double[Ny_co*Nx_co]; // coarse f

	restriction(f_co,res, n_x, n_y,0.125, 0.125 ,0.25 ,1.0/16.0); //full weighted restriction
    delete[] res;
	if(Nx_co==3||Ny_co==3){
		u[0]=f[0]/GS_CENTER; // solve Au=b ??? ToDo: richtig????
	}else{
		double* u_co = new double[Nx_co*Ny_co];
		//for(int k=1;k<nyy;k++)// fuer nyy größer 1
		memset(u_co,0,sizeof(double)*Ny_co*Nx_co);
		initCoarseBD(u, u_co , Nx_co);

		mgm( u_co, f_co, v1, v2,Nx_co ,Ny_co); //recursive call
        delete[] f_co;
		prolongation(u_co,u,n_x,n_y); //prolongation
        delete[] u_co;
	}

	do_gauss_seidel(u,f,n_x,n_y,v2); //post-smoothing*/

}

//calculates residuum
void residuum(double* res,double* f, double* u, const int n_x,const int n_y){
	double hx_local = 1.0/n_x;	
	double hy_local = 1.0/n_y;
	double north, south, west, east = 0;
	for(int j=0;j<n_y;j++){
		for(int i=0;i<n_x;i++){
			west = u[j*n_x+i-1];
			east = u[j*n_x+i+1];	
			north = u[(j+1)*n_x+i];	
			south = u[(j-1)*n_x+i];
			if(i == 0){
				west = 0;	
			}
			if(j == 0){
				south = 0;
			}
			if(i == n_x-1){
				east = 0;
			}
			if(j == n_y-1){
				south = 0;
			}

			res[j*n_x+i] =f[j*n_x+i]- //f-Au
				(1/(hx_local*hy_local))*(GS_HORIZONTAL*(west+ east)+
						GS_VERTICAL*(south+ north)+
						GS_CENTER*u[j*n_x+i]);
		}
	}
}


double calcL2Norm(double *res, int n_x, int n_y){
	double norm = 0;
	for(int j= 0; j<n_y; j++){
		for(int i =0 ; i<n_x; i++){
			norm += res[j*n_y+i]*res[j*n_y+i];
		}
	}
	return sqrt(norm/(n_x*n_y));
}


void measureError(double *u, double gridsize){
	double *error = new double[NX*NY];
	for(int j = 0; j<NY; j++){
		for(int i = 0; i < NX; i++){
			error[j*NX+NY] = sqrt( sin(M_PI*j*H)*sinh(M_PI*i*H)* sin(M_PI*j*H)*sinh(M_PI*i*H)-u[j*NX+i]*u[j*NX+i]);
		}
	}
}
