#include <math.h>
#include <string.h> /*for memcpy function*/
#include "mex.h"
#include "matrix.h"
#include <omp.h>


/*to compile with open mp, use mex residualvel3Dper3_mex.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\LDFLAGS -fopenmp"*/

void transferEtoLdisp(double h, double *x, double *y, double *z, double *X, double *U, int numpoints, double zlength, double xlength, int M, int N, int P, *Edisp, r)
{
    int coef, i, j, k, m,xi, yi, zi, ximin, yimin, zimin;, ximax, yimax, zimax;
    coef =1/(r*r*r);
    double ih, ir;
    ih = 1/h;
    ir = 1/r;

    double xlocal[20][20][20], ylocal[20][20][20], zlocal[20][20][20];
    double xdist[20][20][20], ydist[20][20][20], zdist[20][20][20];


    for(i=0; i<numpoints; i++){
    	xi = mod(ceil(X[3*i]*ih), N);
	yi = ceil(X[3*i+1]*ih);
	zi = mod(ceil(X[3*i+2]*ih), P);

	if(zi == 0 && X[3*i+2]*ih < P-1){
		zi=P-1;}

	if(xi == 0 && X[3*i]*ih < N-1){
		xi = N-1;}
	
	ximin = max(0, xi-10);
	ximax = min(N-2, xi+10);
	yimin = max(0, yi-10);
	yimax = min(M-1, yi+10);
	zimin = max(0, zi-10);
	zimax = min(P-2, zi+10);
	
	for(j=0; j<20; j++){
		for(k=0; k<20; k++){
			for(m=0; m<20; m++){
				xlocal[j][k][m]=h*(ximin+k);
				ylocal[j][k][m]=h*(yimin+j);
				zlocal[j][k][m]=h*(zimin+m);
				xdist[j][k][m]=abs(xlocal[j][k][m]-X[3*i])*ir;
				ydist[j][k][m]=abs(ylocal[j][k][m]-X[3*i+1])*ir;
				zdist[j][k][m]=abs(zlocal[j][k][m]-X[3*i+2])*ir;
				
				}
			}
		}

	xdist = dirac_interp_mex(xdist)*dirac_interp_mex(ydist)*dirac_interp_mex(zdist);
	
	for(j=0; j<20; j++){
		for(k=0; k<20; k++){
			for(m=0; m<20; m++){
				Edisp[yimin+j][ximin+k][zimin+m]=Edisp[yimin+j][ximin+k][zimin+m]+coef*U[i]*xdist[j][k][m];
	
				}
			}
		}
		
		
        }




}
void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
int numpoints;
double h, zlength, xlength, r;
double *x, *y, *z, *X, *U, *Edisp;


/* Input Arguments */
h = mxGetScalar(prhs[0]);
x = mxGetPr(prhs[1]);
y = mxGetPr(prhs[2]);
z = mxGetPr(prhs[3]);
X = mxGetPr(prhs[4]);
U = mxGetPr(prhs[5]);
numpoints = mxGetPr(prhs[6]);
zlength = mxGetPr(prhs[7]);
xlength = mxGetPr(prhs[8]);
r = mxGetScalar(prhs[9])
Edisp = mxGetPr(prhs[10]);

int num = mxGetNumberOfDimensions(prhs[1]);
const int* dims = mxGetDimensions(prhs[1]);

long int M = dims[0];
long int N = dims[1];
long int P = dims[2];

/*start at second element in array (nonboundary)*/
if(omp_get_num_procs()>=4){
	omp_set_num_threads(4);}

void transferEtoLDisp( h, x, y, z, X, U, *numpoints, zlength, xlength, M, N, P, Edisp, r);



plhs[0]=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS, mxREAL);

memcpy(mxGetPr(plhs[0]),Edisp,sizeof(double)*M*N*P);

return;
}
