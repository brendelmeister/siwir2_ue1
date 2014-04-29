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
    NX = NY = getGridPointsDirichlet();
    H = 1.0/(NX-1);


    double* u = new double[NX*NY]; // initialise arrays
    memset(u,0,sizeof(double)*NY*NX);
    initializeGrid(u);
    double* f = new double[NX*NY];
    memset(f,0,sizeof(double)*NY*NX);


    /*do_gauss_seidel(u,f, NX, NY, 100);
      double* res = new double[NY*NX];
      memset(res,0,sizeof(double)*NY*NX);

      residuum(res,f, u, NX,NY);

      double l2norm = calcL2Norm(res, NX, NY);
      cout<<l2norm<<endl;
      save_in_file("residuum_gs100_neue_l2.txt", res, NX, NY);

      memset(res,0,sizeof(double)*NY*NX);
      calculate_L2Norm(res, u, f, NX, NY);
      save_in_file("residuum_gs100_alte_l2.txt", res, NX, NY);
      delete[] res;

     */
    for(int i=0;i<n;i++){ //multigrid steps
        mgm(u,f,2,1,NX, NY);
    }

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

    // LH: Try the Interpolation/Prolongation algorithm from Pflaum-Script p. 6
    int Nx_co=(n_x/2)+1;
    int Ny_co=(n_y/2)+1;

    for(int j = 0; j < Ny_co-1; ++j)
    {
        for(int i = 0; i < Nx_co-1; ++i)
        {
            u_fi[IDX(2*i,2*j)]	    += u_co[j * Nx_co+ i]; // centre 
            u_fi[IDX(2*i+1,2*j)]  	+= 1./2. * (u_co[j * Nx_co+ i] + u_co[j * Nx_co+ i+1]);
            u_fi[IDX(2*i,2*j+1)]    += 1./2. * (u_co[j * Nx_co+ i] + u_co[(j+1) * Nx_co+i]);
            u_fi[IDX(2*i+1,2*j+1)] 	+= 1./4. * (u_co[j * Nx_co+ i] + u_co[(j+1) * Nx_co+i]+ u_co[j * Nx_co+ i+1]+ u_co[(j+1) * Nx_co+ i+1]);
        } 
    }
}


void do_gauss_seidel(double *u, double *f, const int n_x, const int n_y, const int c){

    if(n_x != n_y){
        cout<<"error: grid not quadratic"<<endl;
        exit(EXIT_FAILURE);
    }
    double h = 1.0 / n_x;


    //#pragma omp parallel
    for(int it = 0; it < c; ++it ){

        /*        // gauss seidel "normal" 
                  for(int yi = 1; yi < n_y-1; ++yi){
                  for(int xj = 1; xj < n_x-1; ++xj){
                  u[yi * n_x + xj] = ( h * h * f[yi * n_x + xj]
                  +  u[yi * n_x + xj +1]
                  +  u[yi * n_x + xj -1]
                  +  u[(yi + 1) * n_x + xj]
                  +  u[(yi - 1) * n_x + xj]
                  ) / 4.0;
                  }
                  }
         */


        //red-black
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
    for(int i = 0; i < NX; ++i){
        u[(NY-1) * NX + i] = sin(M_PI*i*H) * sinh(M_PI*(NY-1)*H);
    }
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

// copy boundaries from fine to coarse grid
void initCoarseBD(const double* u_fi, double* u_co, int Nx_co){
    int n_x=Nx_co*2-1;
    for(int i=0;i<Nx_co;i++){
        u_co[0*Nx_co+i] = u_fi[0*2*n_x+i*2];
        u_co[(Nx_co-1)*Nx_co+i] = u_fi[(n_x-1)*n_x+i*2];
    }

    for(int j=0; j<Nx_co;j++){
        u_co[j*Nx_co+0] = u_fi[j*2*n_x+0];
        u_co[j*Nx_co+Nx_co-1] = u_fi[j*2*n_x+n_x-1];
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
    restriction(f_co,res, n_x, n_y); //full weighted restriction

    delete[] res;
    double* c_co=new double[Ny_co*Nx_co]; // coarse f
    memset(c_co,0, sizeof(double)*Ny_co*Nx_co);

    if(Nx_co==3||Ny_co==3)
    {
        //u[0]=f[0]/GS_CENTER; // solve Au=b ??? ToDo: richtig????

        /* LH: Sollte mMn so aussehen (ist im Endeffekt Gauss-Seidel).
         * Drei von den Randwarten sind ja eigentlich 0, lassen wir aber mal trotzdem so stehen, den Code.
         */

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
    do_gauss_seidel(u,f,n_x,n_y,v2); //post-smoothing*/

}

//calculates residuum
void residuum(double* res,double* f, double* u, const int n_x,const int n_y){
    double hx_local = 1.0/n_x;	
    double hy_local = 1.0/n_y;
    for(int j=1;j<n_y-1;j++){
        for(int i=1;i<n_x-1;i++){
            res[IDX(i,j)] =  //f-Au
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

void calculate_L2Norm(double *res, const double *u, const double *f, const int n_x, const int n_y){

    double L2norm = 0.0;
    double sum_res = 0.0;

    double hx_local = 1.0/n_x;	
    double hy_local = 1.0/n_y;
    //    #pragma omp parallel for num_threads(32) if(n_x > 2500 && n_y > 2500)
    for(int yi = 1; yi < n_y-1; ++yi)
    {
        for(int xj = 1; xj < n_x-1; ++xj)
        {
            res[yi * n_x + xj] =  -f[yi * n_x + xj]		
                + (1.0/(hy_local*hx_local)) *(
                        -u[yi * n_x + xj +1]
                        -u[yi * n_x + xj -1]
                        -u[(yi + 1) * n_x + xj]
                        -u[(yi - 1) * n_x + xj]
                        +4.0* u[yi * n_x + xj]);
            sum_res += res[yi * n_x + xj] * res[yi * n_x + xj];
        }
    }
    L2norm = sqrt(sum_res / (n_x  * n_y));
    printf("L2norm: %f\n\n", L2norm);
}

