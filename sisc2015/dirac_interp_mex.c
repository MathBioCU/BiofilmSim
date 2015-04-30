#include <math.h>
#include "mex.h"
#include <string.h>

/*#define r(i,j,k) r[(k+1)*M*N+(j+1)*M+(i+1)]*/

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i, j, k, M, N, P;
    double * r;
   /* M=mxGetNumberOfElements( prhs[0]);*/
   /* int num_dim = mxGetNumberOfDimensions(prhs[0]);*/	
  

  /*Maybe add in error catching for null or empty inputs */

   const int * dims = mxGetDimensions(prhs[0]);
   int num = mxGetNumberOfDimensions(prhs[0]);
/* Catch exceptions if input is not 3D array */
    if(num==3){
    	M = dims[0];
    	N = dims[1];
    	P = dims[2];}
    else if(num==2){
    	M=dims[0];
	N=dims[1];
	P=1;}
    else if(num==1){
    	M=dims[0];
	N=1;
	P=1;}

    r = mxGetPr( prhs[0] );

    for(i=0; i<M*N*P; i++){
            if(r[i]<= 1){
                r[i]=(3-2*r[i]+sqrt(1+4*r[i]-4*r[i]*r[i]))/8;
            }
            else if(r[i]>1 && r[i]<=2){
                r[i]=(5-2*r[i]-sqrt(-7+12*r[i]-4*r[i]*r[i]))/8;
            }
            else{
                r[i]=0;}
	}
    	
   /* for(k=0; k<P-1; k++){
    	for(j=0; j<N-1; j++){
		for(i=0; i<M-1; i++){	
			if(r(i,j,k)<=1){
                		r(i,j,k)=(3-2*r(i,j,k)+sqrt(1+4*r(i,j,k)-4*r(i,j,k)*r(i,j,k)))/8;}
			else if(r(i,j,k)>1 && r(i,j,k)<=2){
				 r(i,j,k)=(5-2*r(i,j,k)-sqrt(-7+12*r(i,j,k)-4*r(i,j,k)*r(i,j,k)))/8;}
			else{ 
				r(i,j,k)=0; 
				}
			}
		}
	}*/
	plhs[0]=mxCreateNumericArray(num,dims,mxDOUBLE_CLASS, mxREAL);
	memcpy(mxGetPr(plhs[0]),r,sizeof(double)*M*N*P);
/*	memcpy(r,r,sizeof(double)*M*N*P); */ /* this line makes it so that there is no output if previous two lines commented out*/
	return;
}
	
