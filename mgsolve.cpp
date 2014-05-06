/****************************************************************************
 *                   FAU Erlangen SS14
 *                   Siwir2 Uebung 1 - Elliptic PDE with Multigrid
 *                   written by Lina Gundelwein, Michael Hildebrand,
 *                   Tim Brendel, Lorenz Hufnagel
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
    NX = NY = pow(2,l)+1;
    H = 1.0/(NX-1);


    double* u = new double[NX*NY]; // initialise arrays
    memset(u,0,sizeof(double)*NY*NX);

    //initializeGrid(u);
    //if neumann
    initBD(u,NX,NY);

    double* f = new double[NX*NY];
    //memset(f,0,sizeof(double)*NY*NX);
    for (int i=0;i<NX*NY;i++)
        f[i]=2.0;

    double* con = new double[NX*NY];
    memset(con,0,sizeof(double)*NY*NX);

    double* resid = new double[NY*NX];
    memset(resid,0,sizeof(double)*NY*NX);

    double l2norm = 0;
    double l2_old = 1.0;

    // time measurements
    struct timeval start, end;
    long seconds, useconds;
    gettimeofday(&start, NULL); 

    for(int i=0;i<n;i++){ 
	//multigrid steps
        mgm(u,f,2,1,NX, NY);

        residuum(resid,f,u,NX,NY);
        // norm and convergence
        l2norm = calcL2Norm(resid, NX,NY);
        cout<<"L2 Norm: "<<l2norm<<endl;

        cout<<"Convergence rate: "<< l2norm / l2_old <<endl;
        l2_old = l2norm;

    }

    gettimeofday(&end, NULL);
    seconds = end.tv_sec -start.tv_sec;
    useconds = end.tv_usec-start.tv_usec;
    if(useconds <0){
        useconds+=1000000;
        seconds--;
    }
    cout<<"Duration: "<<seconds<<" sec "<<useconds<<" usec"<<endl;

    double* error = new double[NX*NY]; 
    save_in_file("solution.txt", u, NX, NY);
    
    measureError(u,error);
    char filename[13];
    sprintf(filename, "error%u.txt", NX-1);
    save_in_file(filename, error , NX, NY);

    delete[] u;
    delete[] f;
    delete[] con;
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

   int Nx_co=(n_x/2)+1;
   int Ny_co=(n_y/2)+1;

   /*
      for(int j = 0; j < Ny_co-1; ++j)
      {
      for(int i = 0; i < Nx_co-1; ++i)
      {
      if (i!=0 && j!=0)
      u_fi[IDX(2*i,2*j)]  	+= u_co[j * Nx_co+ i]; // centre 
      u_fi[IDX(2*i+1,2*j)]	+= 1./2. * (u_co[j * Nx_co+ i] + u_co[j * Nx_co+ i+1]);
      u_fi[IDX(2*i,2*j+1)]    	+= 1./2. * (u_co[j * Nx_co+ i] + u_co[(j+1) * Nx_co+i]);
      u_fi[IDX(2*i+1,2*j+1)] 	+= 1./4. * (u_co[j * Nx_co+ i] + u_co[(j+1) * Nx_co+i]+ u_co[j * Nx_co+ i+1]+ u_co[(j+1) * Nx_co+ i+1]);
      } 
      }
      */

   //set four courners
   int i=0;
   int j=0;
   double center = u_co[j * Nx_co+ i];
   u_fi[IDX(2*i+1,2*j+1)] 	+= 1./4. * center;

   i=Nx_co-1;

   center = u_co[j * Nx_co+ i];
   u_fi[IDX(2*i-1,2*j+1)] 	+= 1./4. * center;

   j=Ny_co-1;
   center = u_co[j * Nx_co+ i];
   u_fi[IDX(2*i-1,2*j-1)] 	+= 1./4. * center;

   i=0;
   center = u_co[j * Nx_co+ i];
   u_fi[IDX(2*i+1,2*j-1)] 	+= 1./4. * center;

   //y-border-columns
   for(j = 1; j < Ny_co-1; ++j)
   {
      i=0;
      center = u_co[j * Nx_co+ i];

      u_fi[IDX(2*i+1,2*j)]	+= 1./2. * center;
      u_fi[IDX(2*i+1,2*j-1)] 	+= 1./4. * center;
      u_fi[IDX(2*i+1,2*j+1)] 	+= 1./4. * center;

      i=Nx_co-1;

      center = u_co[j * Nx_co+ i];

      u_fi[IDX(2*i-1,2*j)]	+= 1./2. * center;
      u_fi[IDX(2*i-1,2*j+1)] 	+= 1./4. * center;
      u_fi[IDX(2*i-1,2*j-1)] 	+= 1./4. * center;
   }

   //x-border-rows
   for(i = 1; i < Nx_co-1; ++i)
   {
      j=0;
      center = u_co[j * Nx_co+ i];


      u_fi[IDX(2*i,2*j+1)]    	+= 1./2. * center;
      u_fi[IDX(2*i-1,2*j+1)] 	+= 1./4. * center;
      u_fi[IDX(2*i+1,2*j+1)] 	+= 1./4. * center;

      j=Ny_co-1;

      center = u_co[j * Nx_co+ i];

      u_fi[IDX(2*i , 2*j-1)]    	+= 1./2. * center;
      u_fi[IDX(2*i+1,2*j-1)] 	+= 1./4. * center;
      u_fi[IDX(2*i-1,2*j-1)] 	+= 1./4. * center;
   }

   //inner grid
   for(j = 1; j < Ny_co-1; ++j)
   {
      for(i = 1; i < Nx_co-1; ++i)
      {

	 center = u_co[j * Nx_co+ i];
	 u_fi[IDX(2*i,2*j)]  	+= center;

	 u_fi[IDX(2*i+1,2*j)]	+= 1./2. * center;
	 u_fi[IDX(2*i-1,2*j)]	+= 1./2. * center;
	 u_fi[IDX(2*i,2*j+1)]    	+= 1./2. * center;
	 u_fi[IDX(2*i,2*j-1)]    	+= 1./2. * center;

	 u_fi[IDX(2*i-1,2*j+1)] 	+= 1./4. * center;
	 u_fi[IDX(2*i+1,2*j-1)] 	+= 1./4. * center;
	 u_fi[IDX(2*i-1,2*j-1)] 	+= 1./4. * center;
	 u_fi[IDX(2*i+1,2*j+1)] 	+= 1./4. * center;
      } 
   }
}


void do_gauss_seidel(double *u, double *f, const int n_x, const int n_y, const int c,const double bdValue){

   if(n_x != n_y){
      cout<<"error: grid not quadratic"<<endl;
      exit(EXIT_FAILURE);
   }
   double h = 1.0 / (n_x-1);

   for(int it = 0; it < c; ++it ){


      //red
      for (int y=1;y<n_y-1;y++)
      {
	 for (int x=(y%2)+1;x<n_x-1;x+=2)
	 {
        u[IDX(x,y)] = 1.0/4.0 * (h*h*f[IDX(x,y)] + (u[IDX(x,y-1)] + u[IDX(x,y+1)] + u[IDX(x-1,y)] + u[IDX(x+1,y)]));
	 }
      }
      //black
      for (int y=1;y<n_y-1;y++)
      {
        for (int x=((y+1)%2)+1;x<n_x-1;x+=2)
        {
            u[IDX(x,y)] = 1.0/4.0 * (h*h*f[IDX(x,y)] + (u[IDX(x,y-1)] + u[IDX(x,y+1)] + u[IDX(x-1,y)] + u[IDX(x+1,y)]));
        }
      }
   }
}

void initializeGrid(double* u){
    for(int i = 0; i < NX; ++i){
        u[(NY-1) * NX + i] = sin(M_PI*i*H) * sinh(M_PI*(NY-1)*H);
    }
}

void initBD(double* u,const int n_x, const int n_y){
   double hx=1./double(n_x-1);

   for(int i = 0; i < n_x; ++i)
   {
      u[IDX(i,0)] = i*hx*(1.0-i*hx);
      u[IDX(i,n_y-1)] = i*hx*(1.0-i*hx);
   }

   setNMBoundary(u,-1.,n_y,n_x);

}


//do restriction from residuum to f_coarse
void restriction(double* f_co,double* res,const int n_x,const int n_y){
    int Nx_co=(n_x/2)+1;
    int Ny_co=(n_y/2)+1;

    //x=0(linker rand) und x=1(rechter rand)
    for(int j=0; j<Ny_co;j++){
        f_co[j*Nx_co+0] = res[j*2*n_x+0];
        f_co[j*Nx_co+Nx_co-1] = res[j*2*n_x+n_x-1];
    }
    //y=0(unterer rand) und y=1(oberer rand)
    for(int i=0; i<Nx_co;i++){
        f_co[0*Nx_co+i] = res[0*2*n_x+i*2];
        f_co[(Ny_co-1)*Nx_co+i] = res[(n_y-1)*n_x+i*2];
    }

    for(int j=1;j<Ny_co-1;j++){
        for(int i=1;i<Nx_co-1;i++){
            f_co[j*Nx_co+i] =
                RES_CENTER*res[IDX(2*i,2*j)]+ //restriction stencil
                RES_HORIZONTAL*(res[(j*2*n_x+i*2)-1]+ res[(j*2*n_x+i*2)+1])+
                RES_VERTICAL*(res[((j*2-1)*n_x+i*2)]+ res[((j*2+1)*n_x+i*2)])+
                RES_CORNER*(res[((j*2-1)*n_x+i*2)-1]+ res[((j*2-1)*n_x+i*2)+1]+
                        res[((j*2+1)*n_x+i*2)-1]+ res[((j*2+1)*n_x+i*2)+1]);
        }
    }
}

// sets NeumanBoundaries at the left and right boundary
void setNMBoundary(double* u,double bdValue,const int n_y,const int n_x){
   double hy = 1./double(n_y);
   for (int y=1;y<n_y-1;y++)
   {
      u[IDX(0,y)] = u[IDX(1,y)]+hy*bdValue;
      u[IDX(n_x-1,y)] = u[IDX(n_x-2,y)]+hy*bdValue;
   }
}


//recursive multigrid function
void mgm(double* u,double* f,int v1,int v2,int n_x, int n_y){


    //if (n_y == NY)
       //setNMBoundary(u,-1.,n_y,n_x);


    do_gauss_seidel(u,f,n_x,n_y,v1,-1.0);//Pre-smoothing

    double* res = new double[n_y*n_x];
    memset(res,0,sizeof(double)*n_y*n_x);


    //if (n_y == NY)
      // setNMBoundary(u,-1.,n_y,n_x);

    residuum(res, f, u, n_x, n_y); //residuum calculation

    int Nx_co=(n_x/2)+1; //calculating coarse grid size
    int Ny_co=(n_y/2)+1;
    double* f_co=new double[Ny_co*Nx_co]; // coarse f
    restriction(f_co,res, n_x, n_y); //full weighted restriction
    delete[] res;
    
    double* c_co=new double[Ny_co*Nx_co]; // coarse f
    memset(c_co,0, sizeof(double)*Ny_co*Nx_co);

    if(Nx_co==3||Ny_co==3)
    {
        double h_co = 1.0/2.;

        c_co[1 * Nx_co + 1] = ( (h_co*h_co) * f[1 * Nx_co + 1]
                +  u[0 * Nx_co + 1]
                +  u[2 * Nx_co + 1]
                +  u[1 * Nx_co + 2]
                +  u[1 * Nx_co + 0]
                ) / 4.0;
    }
    else
    {
        mgm(c_co, f_co, v1, v2,Nx_co ,Ny_co); //recursive call
        delete[] f_co;
    }



    prolongation(c_co,u,n_x,n_y); //prolongation
    delete[] c_co;




     
    do_gauss_seidel(u,f,n_x,n_y,v2,-1.0);//post-smoothing

}

//calculates residuum
void residuum(double* res,double* f, double* u, const int n_x,const int n_y){
    double hx_local = 1.0/(n_x-1);
    double hy_local = 1.0/(n_y-1);


    for(int j=1;j<n_y-1;j++){
        for(int i=1;i<n_x-1;i++){
            res[IDX(i,j)] =  // f-Au
                f[IDX(i,j)] -
                (1.0/(hx_local*hy_local))*
                (4.0*u[IDX(i,j)]
                 - u[IDX( i ,j-1)]
                 - u[IDX( i ,j+1)]
                 - u[IDX(i+1, j )]
                 - u[IDX(i-1, j )])
                ;
        }
    }
}


double calcL2Norm(double *res, int n_x, int n_y){
    double norm = 0.0;
    for(int j= 0; j<n_y; j++){
        for(int i =0 ; i<n_x; i++){
            norm += res[j*n_y+i]*res[j*n_y+i];
        }
    }
    return sqrt(norm/(n_x*n_y));
}


//TODO check what is wrong; change sgn, i think it should be sqrt(u*u-sin*sinh*sin*sinh), but it returns negative
//LH:  i think it should be sqrt(sum((u-sin*sinh)^2))
void measureError(double *u, double *error){
    double h = 1./(NX-1);
    for(int j = 0; j<NY; j++){
        for(int i = 0; i < NX; i++){
            error[j*NX+i] =  u[j*NX+i]-sin(M_PI*i*h)*sinh(M_PI*j*h);
        }
    }
}

