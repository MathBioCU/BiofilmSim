#include <math.h>
#include <string.h> /*for memcpy function*/
#include "mex.h"
#include "matrix.h"
#include <omp.h>


/*to compile with open mp, use mex residualvel3Dper3_mex.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\LDFLAGS -fopenmp"*/

void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
int i,j,k;
double *vhx, *vhy,*vhz,*fhx,*fhy,*fhz, *rx, *ry, *rz;
double *coefx,*coefy,*coefz, *rescoefpm1, *rescoefpm2, *rescoefpm3, *rescoefpp1, *rescoefpp2, *rescoefpp3, *viscpp1, *viscpp2, *viscpp3, *viscpm1, *viscpm2, *viscpm3;
double *rescoefpm1boundF, *rescoefpp1boundF,*rescoefpp2boundF, *rescoefpm2boundF, *rescoefpp3boundF, *rescoefpm3boundF, *rescoefpp1boundR, *rescoefpm1boundR;
double *rescoefpm2boundR, *rescoefpp2boundR, *rescoefpm3boundR, *rescoefpp3boundR, *coefxboundR, *coefyboundR, *coefzboundR, *coefxboundF, *coefyboundF, *coefzboundF;
double *coefxboundFR, *coefyboundFR, *coefzboundFR, *rescoefpp1boundFR, *rescoefpp2boundFR, *rescoefpp3boundFR, *rescoefpm1boundFR, *rescoefpm2boundFR, *rescoefpm3boundFR;
int v1;

/* Input Arguments */
vhx = mxGetPr(prhs[0]);
vhy = mxGetPr(prhs[1]);
vhz = mxGetPr(prhs[2]);
fhx = mxGetPr(prhs[3]);
fhy = mxGetPr(prhs[4]);
fhz = mxGetPr(prhs[5]);
v1 = mxGetScalar(prhs[7]);


coefx = mxGetPr(mxGetField(prhs[6],0,"coefx"));
coefxboundF= mxGetPr(mxGetField(prhs[6],0,"coefxboundF"));
coefxboundR =  mxGetPr(mxGetField(prhs[6],0,"coefxboundR"));
coefxboundFR =  mxGetPr(mxGetField(prhs[6],0,"coefxboundFR"));
coefy =  mxGetPr(mxGetField(prhs[6],0,"coefy"));
coefyboundF =  mxGetPr(mxGetField(prhs[6],0,"coefyboundF"));
coefyboundR =  mxGetPr(mxGetField(prhs[6],0,"coefyboundR"));
coefyboundFR =  mxGetPr(mxGetField(prhs[6],0,"coefyboundFR"));
coefz =  mxGetPr(mxGetField(prhs[6],0,"coefz"));
coefzboundF =  mxGetPr(mxGetField(prhs[6],0,"coefzboundF"));
coefzboundR =  mxGetPr(mxGetField(prhs[6],0,"coefzboundR"));
coefzboundFR =  mxGetPr(mxGetField(prhs[6],0,"coefzboundFR"));
rescoefpp1 =  mxGetPr(mxGetField(prhs[6],0,"rescoefpp1"));
rescoefpp1boundF =  mxGetPr(mxGetField(prhs[6],0,"rescoefpp1boundF"));
rescoefpp1boundR =  mxGetPr(mxGetField(prhs[6],0,"rescoefpp1boundR"));
rescoefpp1boundFR =  mxGetPr(mxGetField(prhs[6],0,"rescoefpp1boundFR"));
rescoefpm1 =  mxGetPr(mxGetField(prhs[6],0,"rescoefpm1"));
rescoefpm1boundF =  mxGetPr(mxGetField(prhs[6],0,"rescoefpm1boundF"));
rescoefpm1boundR = mxGetPr(mxGetField(prhs[6],0,"rescoefpm1boundR"));
rescoefpm1boundFR =  mxGetPr(mxGetField(prhs[6],0,"rescoefpm1boundFR"));
rescoefpp2 =  mxGetPr(mxGetField(prhs[6],0,"rescoefpp2"));
rescoefpp2boundF =  mxGetPr(mxGetField(prhs[6],0,"rescoefpp2boundF"));
rescoefpp2boundR =  mxGetPr(mxGetField(prhs[6],0,"rescoefpp2boundR"));
rescoefpp2boundFR =  mxGetPr(mxGetField(prhs[6],0,"rescoefpp2boundFR"));
rescoefpm2 =  mxGetPr(mxGetField(prhs[6],0,"rescoefpm2"));
rescoefpm2boundF =  mxGetPr(mxGetField(prhs[6],0,"rescoefpm2boundF"));
rescoefpm2boundR =  mxGetPr(mxGetField(prhs[6],0,"rescoefpm2boundR"));
rescoefpm2boundFR =  mxGetPr(mxGetField(prhs[6],0,"rescoefpm2boundFR"));
rescoefpp3 =  mxGetPr(mxGetField(prhs[6],0,"rescoefpp3"));
rescoefpp3boundF =  mxGetPr(mxGetField(prhs[6],0,"rescoefpp3boundF"));
rescoefpp3boundR =  mxGetPr(mxGetField(prhs[6],0,"rescoefpp3boundR"));
rescoefpp3boundFR =  mxGetPr(mxGetField(prhs[6],0,"rescoefpp3boundFR"));
rescoefpm3 =  mxGetPr(mxGetField(prhs[6],0,"rescoefpm3"));
rescoefpm3boundF =  mxGetPr(mxGetField(prhs[6],0,"rescoefpm3boundF"));
rescoefpm3boundR =  mxGetPr(mxGetField(prhs[6],0,"rescoefpm3boundR"));
rescoefpm3boundFR =  mxGetPr(mxGetField(prhs[6],0,"rescoefpm3boundFR"));
viscpp1 =  mxGetPr(mxGetField(prhs[6],0,"viscpp1"));
viscpm1 =  mxGetPr(mxGetField(prhs[6],0,"viscpm1"));
viscpp2 =  mxGetPr(mxGetField(prhs[6],0,"viscpp2"));
viscpm2 =  mxGetPr(mxGetField(prhs[6],0,"viscpm2"));
viscpp3 =  mxGetPr(mxGetField(prhs[6],0,"viscpp3"));
viscpm3  = mxGetPr(mxGetField(prhs[6],0,"viscpm3"));


int num = mxGetNumberOfDimensions(prhs[0]);
const int* dims = mxGetDimensions(prhs[0]);

long int M = dims[0];
long int N = dims[1];
long int P = dims[2];

/*Left Hand Side Outputs
 * Ideally, find a way to act directly on input arrays using pointers rather
 * than reallocating more space for the lhs. This is only because matlab 
 * requires this type of setup.*/



/*start at second element in array (nonboundary)*/

omp_set_num_threads(3);
int p;
for(p=1; p<=v1; p++){
#pragma omp parallel
{
#pragma omp for private(i,j) 
for(k=1; k<=P-2; k=k+2){
	for(j=1; j<=N-2; j=j+2){
		for (i=1; i<=M-2; i=i+2){
			vhx[i+M*j+N*M*k]=(fhx[i+M*j+M*N*k]+2*rescoefpp1[i-1+(M-2)*(j-1)+(k-1)*(M-2)*(N-2)]*vhx[i+(j+1)*M+k*M*N]+
				2*rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+(j-1)*M+k*M*N]+rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+1+j*M+k*N*M]+rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhx[i-1+j*M+k*N*M]+rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+j*M+(k+1)*M*N]+rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+j*M+(k-1)*M*N]
				+viscpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+(j+1)*M+k*N*M]-vhy[i+1+(j-1)*M+k*N*M])+viscpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i-1+(j-1)*M+k*N*M]
				-vhy[i-1+(j+1)*M+k*N*M])+viscpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+(j+1)*M+(k+1)*M*N]-vhz[i+(j-1)*M+(k+1)*M*N])+viscpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhz[i+(j-1)*M+(k-1)*M*N]-vhz[i+(j+1)*M+(k-1)*M*N]))/(coefx[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);
						
			}
		}
	}

#pragma omp for private(i,j)
for(k=1; k<=P-2; k=k+2){
	for(j=2; j<=N-3; j=j+2){
		for (i=2; i<=M-3; i=i+2){
			vhx[i+M*j+N*M*k]=(fhx[i+M*j+M*N*k]+2*rescoefpp1[i-1+(M-2)*(j-1)+(k-1)*(M-2)*(N-2)]*vhx[i+(j+1)*M+k*M*N]+
				2*rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+(j-1)*M+k*M*N]+rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+1+j*M+k*N*M]+rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhx[i-1+j*M+k*N*M]+rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+j*M+(k+1)*M*N]+rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+j*M+(k-1)*M*N]
				+viscpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+(j+1)*M+k*N*M]-vhy[i+1+(j-1)*M+k*N*M])+viscpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i-1+(j-1)*M+k*N*M]
				-vhy[i-1+(j+1)*M+k*N*M])+viscpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+(j+1)*M+(k+1)*M*N]-vhz[i+(j-1)*M+(k+1)*M*N])+viscpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhz[i+(j-1)*M+(k-1)*M*N]-vhz[i+(j+1)*M+(k-1)*M*N]))/(coefx[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);
						
			}
		}
	}


#pragma omp for private(i)
for(k=1; k<=P-2; k=k+2){
	for(i=2; i<=M-3; i=i+2){

		vhx[i+M*(N-1)+N*M*k]=(fhx[i+M*(N-1)+M*N*k]+2*rescoefpp1boundR[i-1+(k-1)*(M-2)]*vhx[i+M+k*M*N]+
			2*rescoefpm1boundR[i-1+(k-1)*(M-2)]*vhx[i+(N-2)*M+k*M*N]+rescoefpp2boundR[i-1+(k-1)*(M-2)]*vhx[i+1+(N-1)*M+k*N*M]+rescoefpm2boundR[i-1+(k-1)*(M-2)]
			*vhx[i-1+(N-1)*M+k*N*M]+rescoefpp3boundR[i-1+(k-1)*(M-2)]*vhx[i+(N-1)*M+(k+1)*M*N]+rescoefpm3boundR[i-1+(k-1)*(M-2)]*vhx[i+(N-1)*M+(k-1)*M*N]
			+viscpp2[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+M+k*N*M]-vhy[i+1+(N-2)*M+k*N*M])+viscpm2[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i-1+(N-2)*M+k*N*M]
			-vhy[i-1+M+k*N*M])+viscpp3[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+M+(k+1)*M*N]-vhz[i+(N-2)*M+(k+1)*M*N])+viscpm3[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]
			*(vhz[i+(N-2)*M+(k-1)*M*N]-vhz[i+M+(k-1)*M*N]))/(coefxboundR[i-1+(k-1)*(M-2)]);
		
		vhx[i+N*M*k]=vhx[i+M*(N-1)+N*M*k];
		}
	}

#pragma omp for private(i,j)
for(k=2; k<=P-3; k=k+2){
	for(j=2; j<=N-3; j=j+2){
		for (i=1; i<=M-2; i=i+2){
			vhx[i+M*j+N*M*k]=(fhx[i+M*j+M*N*k]+2*rescoefpp1[i-1+(M-2)*(j-1)+(k-1)*(M-2)*(N-2)]*vhx[i+(j+1)*M+k*M*N]+
				2*rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+(j-1)*M+k*M*N]+rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+1+j*M+k*N*M]+rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhx[i-1+j*M+k*N*M]+rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+j*M+(k+1)*M*N]+rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+j*M+(k-1)*M*N]
				+viscpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+(j+1)*M+k*N*M]-vhy[i+1+(j-1)*M+k*N*M])+viscpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i-1+(j-1)*M+k*N*M]
				-vhy[i-1+(j+1)*M+k*N*M])+viscpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+(j+1)*M+(k+1)*M*N]-vhz[i+(j-1)*M+(k+1)*M*N])+viscpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhz[i+(j-1)*M+(k-1)*M*N]-vhz[i+(j+1)*M+(k-1)*M*N]))/(coefx[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);
						
			}
		}
	}
#pragma omp for private(i)
for(k=2; k<=P-3; k=k+2){
	for(i=1; i<=M-2; i=i+2){

		vhx[i+M*(N-1)+N*M*k]=(fhx[i+M*(N-1)+M*N*k]+2*rescoefpp1boundR[i-1+(k-1)*(M-2)]*vhx[i+M+k*M*N]+
			2*rescoefpm1boundR[i-1+(k-1)*(M-2)]*vhx[i+(N-2)*M+k*M*N]+rescoefpp2boundR[i-1+(k-1)*(M-2)]*vhx[i+1+(N-1)*M+k*N*M]+rescoefpm2boundR[i-1+(k-1)*(M-2)]
			*vhx[i-1+(N-1)*M+k*N*M]+rescoefpp3boundR[i-1+(k-1)*(M-2)]*vhx[i+(N-1)*M+(k+1)*M*N]+rescoefpm3boundR[i-1+(k-1)*(M-2)]*vhx[i+(N-1)*M+(k-1)*M*N]
			+viscpp2[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+M+k*N*M]-vhy[i+1+(N-2)*M+k*N*M])+viscpm2[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i-1+(N-2)*M+k*N*M]
			-vhy[i-1+M+k*N*M])+viscpp3[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+M+(k+1)*M*N]-vhz[i+(N-2)*M+(k+1)*M*N])+viscpm3[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]
			*(vhz[i+(N-2)*M+(k-1)*M*N]-vhz[i+M+(k-1)*M*N]))/(coefxboundR[i-1+(k-1)*(M-2)]);
		
		vhx[i+N*M*k]=vhx[i+M*(N-1)+N*M*k];
		}
	}
#pragma omp for private(i)
for(j=2; j<=N-3; j=j+2){
	for(i=1; i<=M-2; i=i+2){
		vhx[i+M*j+N*M*(P-1)]=(fhx[i+M*j+M*N*(P-1)]+2*rescoefpp1boundF[i-1+(M-2)*(j-1)]*vhx[i+(j+1)*M+(P-1)*M*N]+
				2*rescoefpm1boundF[i-1+(j-1)*(M-2)]*vhx[i+(j-1)*M+(P-1)*M*N]+rescoefpp2boundF[i-1+(j-1)*(M-2)]*vhx[i+1+j*M+(P-1)*N*M]+rescoefpm2boundF[i-1+(j-1)*(M-2)]
				*vhx[i-1+j*M+(P-1)*N*M]+rescoefpp3boundF[i-1+(j-1)*(M-2)]*vhx[i+j*M+M*N]+rescoefpm3boundF[i-1+(j-1)*(M-2)]*vhx[i+j*M+(P-2)*M*N]
				+viscpp2[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i+1+(j+1)*M+(P-1)*N*M]-vhy[i+1+(j-1)*M+(P-1)*N*M])+viscpm2[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i-1+(j-1)*M+(P-1)*N*M]
				-vhy[i-1+(j+1)*M+(P-1)*N*M])+viscpp3[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhz[i+(j+1)*M+M*N]-vhz[i+(j-1)*M+M*N])+viscpm3[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]
				*(vhz[i+(j-1)*M+(P-2)*M*N]-vhz[i+(j+1)*M+(P-2)*M*N]))/(coefxboundF[i-1+(j-1)*(M-2)]);

		vhx[i+M*j]=vhx[i+M*j+M*N*(P-1)];
		}
	}
#pragma omp for
for(i=1; i<=M-2; i=i+2){             
	vhx[i+M*(N-1)+N*M*(P-1)]=(fhx[i+M*(N-1)+M*N*(P-1)]+2*rescoefpp1boundFR[i-1]*vhx[i+M+(P-1)*M*N]+
		2*rescoefpm1boundFR[i-1]*vhx[i+(N-2)*M+(P-1)*M*N]+rescoefpp2boundFR[i-1]*vhx[i+1+(N-1)*M+(P-1)*N*M]+rescoefpm2boundFR[i-1]
		*vhx[i-1+(N-1)*M+(P-1)*N*M]+rescoefpp3boundFR[i-1]*vhx[i+(N-1)*M+M*N]+rescoefpm3boundFR[i-1]*vhx[i+(N-1)*M+(P-2)*M*N]
		+viscpp2[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i+1+M+(P-1)*N*M]-vhy[i+1+(N-2)*M+(P-1)*N*M])+viscpm2[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i-1+(N-2)*M+(P-1)*N*M]
		-vhy[i-1+M+(P-1)*N*M])+viscpp3[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhz[i+M+M*N]-vhz[i+(N-2)*M+M*N])+viscpm3[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]
		*(vhz[i+(N-2)*M+(P-2)*M*N]-vhz[i+M+(P-2)*M*N]))/(coefxboundFR[i-1]);
	
	vhx[i]=vhx[i+M*(N-1)+N*M*(P-1)];
	vhx[i+M*(N-1)]=vhx[i+M*(N-1)+N*M*(P-1)];
	vhx[i+M*N*(P-1)]=vhx[i+M*(N-1)+N*M*(P-1)];

	}
#pragma omp for private(i,j)
for(k=2; k<=P-3; k=k+2){
	for(j=1; j<=N-2; j=j+2){
		for (i=2; i<=M-3; i=i+2){
			vhx[i+M*j+N*M*k]=(fhx[i+M*j+M*N*k]+2*rescoefpp1[i-1+(M-2)*(j-1)+(k-1)*(M-2)*(N-2)]*vhx[i+(j+1)*M+k*M*N]+
				2*rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+(j-1)*M+k*M*N]+rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+1+j*M+k*N*M]+rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhx[i-1+j*M+k*N*M]+rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+j*M+(k+1)*M*N]+rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+j*M+(k-1)*M*N]
				+viscpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+(j+1)*M+k*N*M]-vhy[i+1+(j-1)*M+k*N*M])+viscpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i-1+(j-1)*M+k*N*M]
				-vhy[i-1+(j+1)*M+k*N*M])+viscpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+(j+1)*M+(k+1)*M*N]-vhz[i+(j-1)*M+(k+1)*M*N])+viscpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhz[i+(j-1)*M+(k-1)*M*N]-vhz[i+(j+1)*M+(k-1)*M*N]))/(coefx[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);
						
			}
		}
	}

#pragma omp for private(i)
for(j=1; j<=N-2; j=j+2){
	for(i=2; i<=M-3; i=i+2){
		vhx[i+M*j+N*M*(P-1)]=(fhx[i+M*j+M*N*(P-1)]+2*rescoefpp1boundF[i-1+(M-2)*(j-1)]*vhx[i+(j+1)*M+(P-1)*M*N]+
				2*rescoefpm1boundF[i-1+(j-1)*(M-2)]*vhx[i+(j-1)*M+(P-1)*M*N]+rescoefpp2boundF[i-1+(j-1)*(M-2)]*vhx[i+1+j*M+(P-1)*N*M]+rescoefpm2boundF[i-1+(j-1)*(M-2)]
				*vhx[i-1+j*M+(P-1)*N*M]+rescoefpp3boundF[i-1+(j-1)*(M-2)]*vhx[i+j*M+M*N]+rescoefpm3boundF[i-1+(j-1)*(M-2)]*vhx[i+j*M+(P-2)*M*N]
				+viscpp2[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i+1+(j+1)*M+(P-1)*N*M]-vhy[i+1+(j-1)*M+(P-1)*N*M])+viscpm2[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i-1+(j-1)*M+(P-1)*N*M]
				-vhy[i-1+(j+1)*M+(P-1)*N*M])+viscpp3[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhz[i+(j+1)*M+M*N]-vhz[i+(j-1)*M+M*N])+viscpm3[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]
				*(vhz[i+(j-1)*M+(P-2)*M*N]-vhz[i+(j+1)*M+(P-2)*M*N]))/(coefxboundF[i-1+(j-1)*(M-2)]);

		vhx[i+M*j]=vhx[i+M*j+M*N*(P-1)];
		}
	}
#pragma omp for private(i,j)
for(k=1; k<=P-2; k=k+2){
	for(j=1; j<=N-2; j=j+2){
		for (i=1; i<=M-2; i=i+2){
			vhy[i+j*M+k*N*M]=(fhy[i+j*M+k*N*M]+rescoefpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+(j+1)*M+k*M*N]+
				rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+(j-1)*M+k*M*N]+2*rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+1+j*M+k*N*M]+2*rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhy[i-1+j*M+k*M*N]+rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+j*M+(k+1)*M*N]+rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+j*M+(k-1)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+1+(j+1)*M+k*N*M]-vhx[i-1+(j+1)*M+k*N*M])+viscpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i-1+(j-1)*M+k*N*M]
				-vhx[i+1+(j-1)*M+k*N*M])+viscpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+1+j*M+(k+1)*M*N]-vhz[i-1+j*M+(k+1)*M*N])+viscpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhz[i-1+j*M+(k-1)*M*N]-vhz[i+1+j*M+(k-1)*M*N]))/(coefy[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);
				}
			}
		}
#pragma omp for private(i,j)
for(k=1; k<=P-2; k=k+2){
	for(j=2; j<=N-3; j=j+2){
		for (i=2; i<=M-3; i=i+2){
			vhy[i+j*M+k*N*M]=(fhy[i+j*M+k*N*M]+rescoefpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+(j+1)*M+k*M*N]+
				rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+(j-1)*M+k*M*N]+2*rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+1+j*M+k*N*M]+2*rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhy[i-1+j*M+k*M*N]+rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+j*M+(k+1)*M*N]+rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+j*M+(k-1)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+1+(j+1)*M+k*N*M]-vhx[i-1+(j+1)*M+k*N*M])+viscpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i-1+(j-1)*M+k*N*M]
				-vhx[i+1+(j-1)*M+k*N*M])+viscpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+1+j*M+(k+1)*M*N]-vhz[i-1+j*M+(k+1)*M*N])+viscpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhz[i-1+j*M+(k-1)*M*N]-vhz[i+1+j*M+(k-1)*M*N]))/(coefy[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);
				}
			}
		}

#pragma omp for private(i)
for(k=1; k<=P-2; k=k+2){
	for(i=2; i<=M-3; i=i+2){
		vhy[i+(N-1)*M+k*N*M]=(fhy[i+(N-1)*M+k*N*M]+rescoefpp1boundR[i-1+(k-1)*(M-2)]*vhy[i+M+k*M*N]+
			rescoefpm1boundR[i-1+(k-1)*(M-2)]*vhy[i+(N-2)*M+k*M*N]+2*rescoefpp2boundR[i-1+(k-1)*(M-2)]*vhy[i+1+(N-1)*M+k*N*M]+2*rescoefpm2boundR[i-1+(k-1)*(M-2)]				
			*vhy[i-1+(N-1)*M+k*M*N]+rescoefpp3boundR[i-1+(k-1)*(M-2)]*vhy[i+(N-1)*M+(k+1)*M*N]+rescoefpm3boundR[i-1+(k-1)*(M-2)]*vhy[i+(N-1)*M+(k-1)*M*N]
			+viscpp1[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+1+M+k*N*M]-vhx[i-1+M+k*N*M])+viscpm1[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i-1+(N-2)*M+k*N*M]
			-vhx[i+1+(N-2)*M+k*N*M])+viscpp3[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+1+(N-1)*M+(k+1)*M*N]-vhz[i-1+(N-1)*M+(k+1)*M*N])+viscpm3[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]
			*(vhz[i-1+(N-1)*M+(k-1)*M*N]-vhz[i+1+(N-1)*M+(k-1)*M*N]))/(coefyboundR[i-1+(k-1)*(M-2)]);
		
		vhy[i+N*M*k]=vhy[i+M*(N-1)+N*M*k];
		}
	}
#pragma omp for private(i,j)
for(k=2; k<=P-3; k=k+2){
	for(j=2; j<=N-3; j=j+2){
		for (i=1; i<=M-2; i=i+2){
			vhy[i+j*M+k*N*M]=(fhy[i+j*M+k*N*M]+rescoefpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+(j+1)*M+k*M*N]+
				rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+(j-1)*M+k*M*N]+2*rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+1+j*M+k*N*M]+2*rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhy[i-1+j*M+k*M*N]+rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+j*M+(k+1)*M*N]+rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+j*M+(k-1)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+1+(j+1)*M+k*N*M]-vhx[i-1+(j+1)*M+k*N*M])+viscpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i-1+(j-1)*M+k*N*M]
				-vhx[i+1+(j-1)*M+k*N*M])+viscpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+1+j*M+(k+1)*M*N]-vhz[i-1+j*M+(k+1)*M*N])+viscpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhz[i-1+j*M+(k-1)*M*N]-vhz[i+1+j*M+(k-1)*M*N]))/(coefy[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);
				}
			}
		}
#pragma omp for private(i)
for(k=2; k<=P-3; k=k+2){
	for(i=1; i<=M-2; i=i+2){
		vhy[i+(N-1)*M+k*N*M]=(fhy[i+(N-1)*M+k*N*M]+rescoefpp1boundR[i-1+(k-1)*(M-2)]*vhy[i+M+k*M*N]+
			rescoefpm1boundR[i-1+(k-1)*(M-2)]*vhy[i+(N-2)*M+k*M*N]+2*rescoefpp2boundR[i-1+(k-1)*(M-2)]*vhy[i+1+(N-1)*M+k*N*M]+2*rescoefpm2boundR[i-1+(k-1)*(M-2)]				
			*vhy[i-1+(N-1)*M+k*M*N]+rescoefpp3boundR[i-1+(k-1)*(M-2)]*vhy[i+(N-1)*M+(k+1)*M*N]+rescoefpm3boundR[i-1+(k-1)*(M-2)]*vhy[i+(N-1)*M+(k-1)*M*N]
			+viscpp1[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+1+M+k*N*M]-vhx[i-1+M+k*N*M])+viscpm1[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i-1+(N-2)*M+k*N*M]
			-vhx[i+1+(N-2)*M+k*N*M])+viscpp3[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+1+(N-1)*M+(k+1)*M*N]-vhz[i-1+(N-1)*M+(k+1)*M*N])+viscpm3[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]
			*(vhz[i-1+(N-1)*M+(k-1)*M*N]-vhz[i+1+(N-1)*M+(k-1)*M*N]))/(coefyboundR[i-1+(k-1)*(M-2)]);
		
		vhy[i+N*M*k]=vhy[i+M*(N-1)+N*M*k];
		}
	}
#pragma omp for private(i)
for(j=2; j<=N-3; j=j+2){
	for(i=1; i<=M-2; i=i+2){
		vhy[i+j*M+(P-1)*N*M]=(fhy[i+j*M+(P-1)*N*M]+rescoefpp1boundF[i-1+(j-1)*(M-2)]*vhy[i+(j+1)*M+(P-1)*M*N]+
				rescoefpm1boundF[i-1+(j-1)*(M-2)]*vhy[i+(j-1)*M+(P-1)*M*N]+2*rescoefpp2boundF[i-1+(j-1)*(M-2)]*vhy[i+1+j*M+(P-1)*N*M]+2*rescoefpm2boundF[i-1+(j-1)*(M-2)]
				*vhy[i-1+j*M+(P-1)*M*N]+rescoefpp3boundF[i-1+(j-1)*(M-2)]*vhy[i+j*M+M*N]+rescoefpm3boundF[i-1+(j-1)*(M-2)]*vhy[i+j*M+(P-2)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+1+(j+1)*M+(P-1)*N*M]-vhx[i-1+(j+1)*M+(P-1)*N*M])+viscpm1[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i-1+(j-1)*M+(P-1)*N*M]
				-vhx[i+1+(j-1)*M+(P-1)*N*M])+viscpp3[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhz[i+1+j*M+M*N]-vhz[i-1+j*M+M*N])+viscpm3[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]
				*(vhz[i-1+j*M+(P-2)*M*N]-vhz[i+1+j*M+(P-2)*M*N]))/(coefyboundF[i-1+(j-1)*(M-2)]);
				
		vhy[i+M*j]=vhy[i+M*j+M*N*(P-1)];
		}
	}

#pragma omp for
for(i=1; i<=M-2; i=i+2){
	vhy[i+(N-1)*M+(P-1)*N*M]=(fhy[i+(N-1)*M+(P-1)*N*M]+rescoefpp1boundFR[i-1]*vhy[i+M+(P-1)*M*N]+
		rescoefpm1boundFR[i-1]*vhy[i+(N-2)*M+(P-1)*M*N]+2*rescoefpp2boundFR[i-1]*vhy[i+1+(N-1)*M+(P-1)*N*M]+2*rescoefpm2boundFR[i-1]				
		*vhy[i-1+(N-1)*M+(P-1)*M*N]+rescoefpp3boundFR[i-1]*vhy[i+(N-1)*M+M*N]+rescoefpm3boundFR[i-1]*vhy[i+(N-1)*M+(P-2)*M*N]
		+viscpp1[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+1+M+(P-1)*N*M]-vhx[i-1+M+(P-1)*N*M])+viscpm1[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i-1+(N-2)*M+(P-1)*N*M]
		-vhx[i+1+(N-2)*M+(P-1)*N*M])+viscpp3[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhz[i+1+(N-1)*M+M*N]-vhz[i-1+(N-1)*M+M*N])+viscpm3[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]
		*(vhz[i-1+(N-1)*M+(P-2)*M*N]-vhz[i+1+(N-1)*M+(P-2)*M*N]))/(coefyboundFR[i-1]);

	vhy[i]=vhy[i+M*(N-1)+N*M*(P-1)];
	vhy[i+M*(N-1)]=vhy[i+M*(N-1)+N*M*(P-1)];
	vhy[i+M*N*(P-1)]=vhy[i+M*(N-1)+N*M*(P-1)];
	}
#pragma omp for private(i,j)
for(k=2; k<=P-3; k=k+2){
	for(j=1; j<=N-2; j=j+2){
		for (i=2; i<=M-3; i=i+2){
			vhy[i+j*M+k*N*M]=(fhy[i+j*M+k*N*M]+rescoefpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+(j+1)*M+k*M*N]+
				rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+(j-1)*M+k*M*N]+2*rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+1+j*M+k*N*M]+2*rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhy[i-1+j*M+k*M*N]+rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+j*M+(k+1)*M*N]+rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+j*M+(k-1)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+1+(j+1)*M+k*N*M]-vhx[i-1+(j+1)*M+k*N*M])+viscpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i-1+(j-1)*M+k*N*M]
				-vhx[i+1+(j-1)*M+k*N*M])+viscpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+1+j*M+(k+1)*M*N]-vhz[i-1+j*M+(k+1)*M*N])+viscpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhz[i-1+j*M+(k-1)*M*N]-vhz[i+1+j*M+(k-1)*M*N]))/(coefy[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);
				}
			}
		}


#pragma omp for private(i)
for(j=1; j<=N-2; j=j+2){
	for(i=2; i<=M-3; i=i+2){
		vhy[i+j*M+(P-1)*N*M]=(fhy[i+j*M+(P-1)*N*M]+rescoefpp1boundF[i-1+(j-1)*(M-2)]*vhy[i+(j+1)*M+(P-1)*M*N]+
				rescoefpm1boundF[i-1+(j-1)*(M-2)]*vhy[i+(j-1)*M+(P-1)*M*N]+2*rescoefpp2boundF[i-1+(j-1)*(M-2)]*vhy[i+1+j*M+(P-1)*N*M]+2*rescoefpm2boundF[i-1+(j-1)*(M-2)]
				*vhy[i-1+j*M+(P-1)*M*N]+rescoefpp3boundF[i-1+(j-1)*(M-2)]*vhy[i+j*M+M*N]+rescoefpm3boundF[i-1+(j-1)*(M-2)]*vhy[i+j*M+(P-2)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+1+(j+1)*M+(P-1)*N*M]-vhx[i-1+(j+1)*M+(P-1)*N*M])+viscpm1[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i-1+(j-1)*M+(P-1)*N*M]
				-vhx[i+1+(j-1)*M+(P-1)*N*M])+viscpp3[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhz[i+1+j*M+M*N]-vhz[i-1+j*M+M*N])+viscpm3[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]
				*(vhz[i-1+j*M+(P-2)*M*N]-vhz[i+1+j*M+(P-2)*M*N]))/(coefyboundF[i-1+(j-1)*(M-2)]);
				
		vhy[i+M*j]=vhy[i+M*j+M*N*(P-1)];
		}
	}
#pragma omp for private(i,j)
for(k=1; k<=P-2; k=k+2){
	for(j=1; j<=N-2; j=j+2){
		for(i=1; i<=M-2; i=i+2){
			vhz[i+j*M+k*N*M]=(fhz[i+j*M+k*N*M]+rescoefpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+(j+1)*M+k*M*N]+
				rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+(j-1)*M+k*M*N]+rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+1+j*M+k*N*M]+rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhz[i-1+j*M+k*M*N]+2*rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+j*M+(k+1)*M*N]+2*rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+j*M+(k-1)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(j+1)*M+(k+1)*N*M]-vhx[i+(j+1)*M+(k-1)*N*M])+viscpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(j-1)*M+(k-1)*N*M]
				-vhx[i+(j-1)*M+(k+1)*N*M])+viscpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+j*M+(k+1)*M*N]-vhy[i+1+j*M+(k-1)*M*N])+viscpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhy[i-1+j*M+(k-1)*M*N]-vhy[i-1+j*M+(k+1)*M*N]))/(coefz[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);

			}
		}
	}

#pragma omp for private(i,j)
for(k=1; k<=P-2; k=k+2){
	for(j=2; j<=N-3; j=j+2){
		for(i=2; i<=M-3; i=i+2){
			vhz[i+j*M+k*N*M]=(fhz[i+j*M+k*N*M]+rescoefpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+(j+1)*M+k*M*N]+
				rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+(j-1)*M+k*M*N]+rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+1+j*M+k*N*M]+rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhz[i-1+j*M+k*M*N]+2*rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+j*M+(k+1)*M*N]+2*rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+j*M+(k-1)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(j+1)*M+(k+1)*N*M]-vhx[i+(j+1)*M+(k-1)*N*M])+viscpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(j-1)*M+(k-1)*N*M]
				-vhx[i+(j-1)*M+(k+1)*N*M])+viscpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+j*M+(k+1)*M*N]-vhy[i+1+j*M+(k-1)*M*N])+viscpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhy[i-1+j*M+(k-1)*M*N]-vhy[i-1+j*M+(k+1)*M*N]))/(coefz[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);

			}
		}
	}
#pragma omp for private(i)
for(k=1; k<=P-2; k=k+2){
	for(i=2; i<=M-3; i=i+2){
		vhz[i+(N-1)*M+k*N*M]=(fhz[i+(N-1)*M+k*N*M]+rescoefpp1boundR[i-1+(k-1)*(M-2)]*vhz[i+M+k*M*N]+
			rescoefpm1boundR[i-1+(k-1)*(M-2)]*vhz[i+(N-2)*M+k*M*N]+rescoefpp2boundR[i-1+(k-1)*(M-2)]*vhz[i+1+(N-1)*M+k*N*M]+rescoefpm2boundR[i-1+(k-1)*(M-2)]
			*vhz[i-1+(N-1)*M+k*M*N]+2*rescoefpp3boundR[i-1+(k-1)*(M-2)]*vhz[i+(N-1)*M+(k+1)*M*N]+2*rescoefpm3boundR[i-1+(k-1)*(M-2)]*vhz[i+(N-1)*M+(k-1)*M*N]
			+viscpp1[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+M+(k+1)*N*M]-vhx[i+M+(k-1)*N*M])+viscpm1[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(N-2)*M+(k-1)*N*M]
			-vhx[i+(N-2)*M+(k+1)*N*M])+viscpp2[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+(N-1)*M+(k+1)*M*N]-vhy[i+1+(N-1)*M+(k-1)*M*N])+viscpm2[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]
			*(vhy[i-1+(N-1)*M+(k-1)*M*N]-vhy[i-1+(N-1)*M+(k+1)*M*N]))/(coefzboundR[i-1+(k-1)*(M-2)]);

		vhz[i+N*M*k]=vhz[i+M*(N-1)+N*M*k];
		}
	}
#pragma omp for private(i,j)
for(k=2; k<=P-3; k=k+2){
	for(j=2; j<=N-3; j=j+2){
		for(i=1; i<=M-2; i=i+2){
			vhz[i+j*M+k*N*M]=(fhz[i+j*M+k*N*M]+rescoefpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+(j+1)*M+k*M*N]+
				rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+(j-1)*M+k*M*N]+rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+1+j*M+k*N*M]+rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhz[i-1+j*M+k*M*N]+2*rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+j*M+(k+1)*M*N]+2*rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+j*M+(k-1)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(j+1)*M+(k+1)*N*M]-vhx[i+(j+1)*M+(k-1)*N*M])+viscpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(j-1)*M+(k-1)*N*M]
				-vhx[i+(j-1)*M+(k+1)*N*M])+viscpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+j*M+(k+1)*M*N]-vhy[i+1+j*M+(k-1)*M*N])+viscpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhy[i-1+j*M+(k-1)*M*N]-vhy[i-1+j*M+(k+1)*M*N]))/(coefz[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);

			}
		}
	}
#pragma omp for private(i)
for(k=2; k<=P-3; k=k+2){
	for(i=1; i<=M-2; i=i+2){
		vhz[i+(N-1)*M+k*N*M]=(fhz[i+(N-1)*M+k*N*M]+rescoefpp1boundR[i-1+(k-1)*(M-2)]*vhz[i+M+k*M*N]+
			rescoefpm1boundR[i-1+(k-1)*(M-2)]*vhz[i+(N-2)*M+k*M*N]+rescoefpp2boundR[i-1+(k-1)*(M-2)]*vhz[i+1+(N-1)*M+k*N*M]+rescoefpm2boundR[i-1+(k-1)*(M-2)]
			*vhz[i-1+(N-1)*M+k*M*N]+2*rescoefpp3boundR[i-1+(k-1)*(M-2)]*vhz[i+(N-1)*M+(k+1)*M*N]+2*rescoefpm3boundR[i-1+(k-1)*(M-2)]*vhz[i+(N-1)*M+(k-1)*M*N]
			+viscpp1[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+M+(k+1)*N*M]-vhx[i+M+(k-1)*N*M])+viscpm1[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(N-2)*M+(k-1)*N*M]
			-vhx[i+(N-2)*M+(k+1)*N*M])+viscpp2[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+(N-1)*M+(k+1)*M*N]-vhy[i+1+(N-1)*M+(k-1)*M*N])+viscpm2[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]
			*(vhy[i-1+(N-1)*M+(k-1)*M*N]-vhy[i-1+(N-1)*M+(k+1)*M*N]))/(coefzboundR[i-1+(k-1)*(M-2)]);

		vhz[i+N*M*k]=vhz[i+M*(N-1)+N*M*k];
		}
	}
#pragma omp for private(i)
for(j=2; j<=N-3; j=j+2){
	for(i=1; i<=M-2; i=i+2){
		vhz[i+j*M+(P-1)*N*M]=(fhz[i+j*M+(P-1)*N*M]+rescoefpp1boundF[i-1+(j-1)*(M-2)]*vhz[i+(j+1)*M+(P-1)*M*N]+
				rescoefpm1boundF[i-1+(j-1)*(M-2)]*vhz[i+(j-1)*M+(P-1)*M*N]+rescoefpp2boundF[i-1+(j-1)*(M-2)]*vhz[i+1+j*M+(P-1)*N*M]+rescoefpm2boundF[i-1+(j-1)*(M-2)]
				*vhz[i-1+j*M+(P-1)*M*N]+2*rescoefpp3boundF[i-1+(j-1)*(M-2)]*vhz[i+j*M+M*N]+2*rescoefpm3boundF[i-1+(j-1)*(M-2)]*vhz[i+j*M+M*N*(P-2)]
				+viscpp1[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+(j+1)*M+M*N]-vhx[i+(j+1)*M+(P-2)*N*M])+viscpm1[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+(j-1)*M+(P-2)*N*M]
				-vhx[i+(j-1)*M+M*N])+viscpp2[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i+1+j*M+M*N]-vhy[i+1+j*M+(P-2)*M*N])+viscpm2[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]
				*(vhy[i-1+j*M+(P-2)*M*N]-vhy[i-1+j*M+M*N]))/(coefzboundF[i-1+(j-1)*(M-2)]);

		vhz[i+M*j]=vhz[i+M*j+M*N*(P-1)];
		}
	}

#pragma omp for
for(i=1; i<=M-2; i=i+2){
	vhz[i+(N-1)*M+(P-1)*N*M]=(fhz[i+(N-1)*M+(P-1)*N*M]+rescoefpp1boundFR[i-1]*vhz[i+M+(P-1)*M*N]+
		rescoefpm1boundFR[i-1]*vhz[i+(N-2)*M+(P-1)*M*N]+rescoefpp2boundFR[i-1]*vhz[i+1+(N-1)*M+(P-1)*N*M]+rescoefpm2boundFR[i-1]
		*vhz[i-1+(N-1)*M+(P-1)*M*N]+2*rescoefpp3boundFR[i-1]*vhz[i+(N-1)*M+M*N]+2*rescoefpm3boundFR[i-1]*vhz[i+(N-1)*M+(P-2)*M*N]
		+viscpp1[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+M+M*N]-vhx[i+M+(P-2)*N*M])+viscpm1[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+(N-2)*M+(P-2)*N*M]
		-vhx[i+(N-2)*M+M*N])+viscpp2[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i+1+(N-1)*M+M*N]-vhy[i+1+(N-1)*M+(P-2)*M*N])+viscpm2[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]
		*(vhy[i-1+(N-1)*M+(P-2)*M*N]-vhy[i-1+(N-1)*M+M*N]))/(coefzboundFR[i-1]);

	vhz[i]=vhz[i+M*(N-1)+N*M*(P-1)];
	vhz[i+M*(N-1)]=vhz[i+M*(N-1)+N*M*(P-1)];
	vhz[i+M*N*(P-1)]=vhz[i+M*(N-1)+N*M*(P-1)];

	}
#pragma omp for private(i,j)
for(k=2; k<=P-3; k=k+2){
	for(j=1; j<=N-2; j=j+2){
		for(i=2; i<=M-3; i=i+2){
			vhz[i+j*M+k*N*M]=(fhz[i+j*M+k*N*M]+rescoefpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+(j+1)*M+k*M*N]+
				rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+(j-1)*M+k*M*N]+rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+1+j*M+k*N*M]+rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhz[i-1+j*M+k*M*N]+2*rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+j*M+(k+1)*M*N]+2*rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+j*M+(k-1)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(j+1)*M+(k+1)*N*M]-vhx[i+(j+1)*M+(k-1)*N*M])+viscpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(j-1)*M+(k-1)*N*M]
				-vhx[i+(j-1)*M+(k+1)*N*M])+viscpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+j*M+(k+1)*M*N]-vhy[i+1+j*M+(k-1)*M*N])+viscpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhy[i-1+j*M+(k-1)*M*N]-vhy[i-1+j*M+(k+1)*M*N]))/(coefz[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);

			}
		}
	}


#pragma omp for private(i)
for(j=1; j<=N-2; j=j+2){
	for(i=2; i<=M-3; i=i+2){
		vhz[i+j*M+(P-1)*N*M]=(fhz[i+j*M+(P-1)*N*M]+rescoefpp1boundF[i-1+(j-1)*(M-2)]*vhz[i+(j+1)*M+(P-1)*M*N]+
				rescoefpm1boundF[i-1+(j-1)*(M-2)]*vhz[i+(j-1)*M+(P-1)*M*N]+rescoefpp2boundF[i-1+(j-1)*(M-2)]*vhz[i+1+j*M+(P-1)*N*M]+rescoefpm2boundF[i-1+(j-1)*(M-2)]
				*vhz[i-1+j*M+(P-1)*M*N]+2*rescoefpp3boundF[i-1+(j-1)*(M-2)]*vhz[i+j*M+M*N]+2*rescoefpm3boundF[i-1+(j-1)*(M-2)]*vhz[i+j*M+M*N*(P-2)]
				+viscpp1[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+(j+1)*M+M*N]-vhx[i+(j+1)*M+(P-2)*N*M])+viscpm1[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+(j-1)*M+(P-2)*N*M]
				-vhx[i+(j-1)*M+M*N])+viscpp2[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i+1+j*M+M*N]-vhy[i+1+j*M+(P-2)*M*N])+viscpm2[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]
				*(vhy[i-1+j*M+(P-2)*M*N]-vhy[i-1+j*M+M*N]))/(coefzboundF[i-1+(j-1)*(M-2)]);

		vhz[i+M*j]=vhz[i+M*j+M*N*(P-1)];
		}
	}
#pragma omp for private(i,j)
for(k=1; k<=P-2; k=k+2){
	for(j=2; j<=N-3; j=j+2){
		for (i=1; i<=M-2; i=i+2){
			vhx[i+M*j+N*M*k]=(fhx[i+M*j+M*N*k]+2*rescoefpp1[i-1+(M-2)*(j-1)+(k-1)*(M-2)*(N-2)]*vhx[i+(j+1)*M+k*M*N]+
				2*rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+(j-1)*M+k*M*N]+rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+1+j*M+k*N*M]+rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhx[i-1+j*M+k*N*M]+rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+j*M+(k+1)*M*N]+rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+j*M+(k-1)*M*N]
				+viscpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+(j+1)*M+k*N*M]-vhy[i+1+(j-1)*M+k*N*M])+viscpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i-1+(j-1)*M+k*N*M]
				-vhy[i-1+(j+1)*M+k*N*M])+viscpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+(j+1)*M+(k+1)*M*N]-vhz[i+(j-1)*M+(k+1)*M*N])+viscpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhz[i+(j-1)*M+(k-1)*M*N]-vhz[i+(j+1)*M+(k-1)*M*N]))/(coefx[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);
						
			}
		}
	}
#pragma omp for private(i)
for(k=1; k<=P-2; k=k+2){
	for(i=1; i<=M-2; i=i+2){

		vhx[i+M*(N-1)+N*M*k]=(fhx[i+M*(N-1)+M*N*k]+2*rescoefpp1boundR[i-1+(k-1)*(M-2)]*vhx[i+M+k*M*N]+
			2*rescoefpm1boundR[i-1+(k-1)*(M-2)]*vhx[i+(N-2)*M+k*M*N]+rescoefpp2boundR[i-1+(k-1)*(M-2)]*vhx[i+1+(N-1)*M+k*N*M]+rescoefpm2boundR[i-1+(k-1)*(M-2)]
			*vhx[i-1+(N-1)*M+k*N*M]+rescoefpp3boundR[i-1+(k-1)*(M-2)]*vhx[i+(N-1)*M+(k+1)*M*N]+rescoefpm3boundR[i-1+(k-1)*(M-2)]*vhx[i+(N-1)*M+(k-1)*M*N]
			+viscpp2[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+M+k*N*M]-vhy[i+1+(N-2)*M+k*N*M])+viscpm2[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i-1+(N-2)*M+k*N*M]
			-vhy[i-1+M+k*N*M])+viscpp3[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+M+(k+1)*M*N]-vhz[i+(N-2)*M+(k+1)*M*N])+viscpm3[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]
			*(vhz[i+(N-2)*M+(k-1)*M*N]-vhz[i+M+(k-1)*M*N]))/(coefxboundR[i-1+(k-1)*(M-2)]);
		
		vhx[i+N*M*k]=vhx[i+M*(N-1)+N*M*k];
		}
	}
#pragma omp for private(i,j)
for(k=1; k<=P-2; k=k+2){
	for(j=1; j<=N-2; j=j+2){
		for (i=2; i<=M-3; i=i+2){
			vhx[i+M*j+N*M*k]=(fhx[i+M*j+M*N*k]+2*rescoefpp1[i-1+(M-2)*(j-1)+(k-1)*(M-2)*(N-2)]*vhx[i+(j+1)*M+k*M*N]+
				2*rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+(j-1)*M+k*M*N]+rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+1+j*M+k*N*M]+rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhx[i-1+j*M+k*N*M]+rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+j*M+(k+1)*M*N]+rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+j*M+(k-1)*M*N]
				+viscpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+(j+1)*M+k*N*M]-vhy[i+1+(j-1)*M+k*N*M])+viscpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i-1+(j-1)*M+k*N*M]
				-vhy[i-1+(j+1)*M+k*N*M])+viscpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+(j+1)*M+(k+1)*M*N]-vhz[i+(j-1)*M+(k+1)*M*N])+viscpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhz[i+(j-1)*M+(k-1)*M*N]-vhz[i+(j+1)*M+(k-1)*M*N]))/(coefx[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);
						
			}
		}
	}
#pragma omp for private(i,j)
for(k=2; k<=P-3; k=k+2){
	for(j=1; j<=N-2; j=j+2){
		for (i=1; i<=M-2; i=i+2){
			vhx[i+M*j+N*M*k]=(fhx[i+M*j+M*N*k]+2*rescoefpp1[i-1+(M-2)*(j-1)+(k-1)*(M-2)*(N-2)]*vhx[i+(j+1)*M+k*M*N]+
				2*rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+(j-1)*M+k*M*N]+rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+1+j*M+k*N*M]+rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhx[i-1+j*M+k*N*M]+rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+j*M+(k+1)*M*N]+rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+j*M+(k-1)*M*N]
				+viscpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+(j+1)*M+k*N*M]-vhy[i+1+(j-1)*M+k*N*M])+viscpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i-1+(j-1)*M+k*N*M]
				-vhy[i-1+(j+1)*M+k*N*M])+viscpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+(j+1)*M+(k+1)*M*N]-vhz[i+(j-1)*M+(k+1)*M*N])+viscpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhz[i+(j-1)*M+(k-1)*M*N]-vhz[i+(j+1)*M+(k-1)*M*N]))/(coefx[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);
						
			}
		}
	}
#pragma omp for private(i)
for(j=1; j<=N-2; j=j+2){
	for(i=1; i<=M-2; i=i+2){
		vhx[i+M*j+N*M*(P-1)]=(fhx[i+M*j+M*N*(P-1)]+2*rescoefpp1boundF[i-1+(M-2)*(j-1)]*vhx[i+(j+1)*M+(P-1)*M*N]+
				2*rescoefpm1boundF[i-1+(j-1)*(M-2)]*vhx[i+(j-1)*M+(P-1)*M*N]+rescoefpp2boundF[i-1+(j-1)*(M-2)]*vhx[i+1+j*M+(P-1)*N*M]+rescoefpm2boundF[i-1+(j-1)*(M-2)]
				*vhx[i-1+j*M+(P-1)*N*M]+rescoefpp3boundF[i-1+(j-1)*(M-2)]*vhx[i+j*M+M*N]+rescoefpm3boundF[i-1+(j-1)*(M-2)]*vhx[i+j*M+(P-2)*M*N]
				+viscpp2[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i+1+(j+1)*M+(P-1)*N*M]-vhy[i+1+(j-1)*M+(P-1)*N*M])+viscpm2[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i-1+(j-1)*M+(P-1)*N*M]
				-vhy[i-1+(j+1)*M+(P-1)*N*M])+viscpp3[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhz[i+(j+1)*M+M*N]-vhz[i+(j-1)*M+M*N])+viscpm3[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]
				*(vhz[i+(j-1)*M+(P-2)*M*N]-vhz[i+(j+1)*M+(P-2)*M*N]))/(coefxboundF[i-1+(j-1)*(M-2)]);

		vhx[i+M*j]=vhx[i+M*j+M*N*(P-1)];
		}
	}
#pragma omp for private(i,j)
for(k=2; k<=P-3; k=k+2){
	for(j=2; j<=N-3; j=j+2){
		for (i=2; i<=M-3; i=i+2){
			vhx[i+M*j+N*M*k]=(fhx[i+M*j+M*N*k]+2*rescoefpp1[i-1+(M-2)*(j-1)+(k-1)*(M-2)*(N-2)]*vhx[i+(j+1)*M+k*M*N]+
				2*rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+(j-1)*M+k*M*N]+rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+1+j*M+k*N*M]+rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhx[i-1+j*M+k*N*M]+rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+j*M+(k+1)*M*N]+rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+j*M+(k-1)*M*N]
				+viscpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+(j+1)*M+k*N*M]-vhy[i+1+(j-1)*M+k*N*M])+viscpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i-1+(j-1)*M+k*N*M]
				-vhy[i-1+(j+1)*M+k*N*M])+viscpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+(j+1)*M+(k+1)*M*N]-vhz[i+(j-1)*M+(k+1)*M*N])+viscpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhz[i+(j-1)*M+(k-1)*M*N]-vhz[i+(j+1)*M+(k-1)*M*N]))/(coefx[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);
						
			}
		}
	}


#pragma omp for private(i)
for(k=2; k<=P-3; k=k+2){
	for(i=2; i<=M-3; i=i+2){

		vhx[i+M*(N-1)+N*M*k]=(fhx[i+M*(N-1)+M*N*k]+2*rescoefpp1boundR[i-1+(k-1)*(M-2)]*vhx[i+M+k*M*N]+
			2*rescoefpm1boundR[i-1+(k-1)*(M-2)]*vhx[i+(N-2)*M+k*M*N]+rescoefpp2boundR[i-1+(k-1)*(M-2)]*vhx[i+1+(N-1)*M+k*N*M]+rescoefpm2boundR[i-1+(k-1)*(M-2)]
			*vhx[i-1+(N-1)*M+k*N*M]+rescoefpp3boundR[i-1+(k-1)*(M-2)]*vhx[i+(N-1)*M+(k+1)*M*N]+rescoefpm3boundR[i-1+(k-1)*(M-2)]*vhx[i+(N-1)*M+(k-1)*M*N]
			+viscpp2[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+M+k*N*M]-vhy[i+1+(N-2)*M+k*N*M])+viscpm2[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i-1+(N-2)*M+k*N*M]
			-vhy[i-1+M+k*N*M])+viscpp3[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+M+(k+1)*M*N]-vhz[i+(N-2)*M+(k+1)*M*N])+viscpm3[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]
			*(vhz[i+(N-2)*M+(k-1)*M*N]-vhz[i+M+(k-1)*M*N]))/(coefxboundR[i-1+(k-1)*(M-2)]);
		
		vhx[i+N*M*k]=vhx[i+M*(N-1)+N*M*k];
		}
	}
#pragma omp for private(i)
for(j=2; j<=N-3; j=j+2){
	for(i=2; i<=M-3; i=i+2){
		vhx[i+M*j+N*M*(P-1)]=(fhx[i+M*j+M*N*(P-1)]+2*rescoefpp1boundF[i-1+(M-2)*(j-1)]*vhx[i+(j+1)*M+(P-1)*M*N]+
				2*rescoefpm1boundF[i-1+(j-1)*(M-2)]*vhx[i+(j-1)*M+(P-1)*M*N]+rescoefpp2boundF[i-1+(j-1)*(M-2)]*vhx[i+1+j*M+(P-1)*N*M]+rescoefpm2boundF[i-1+(j-1)*(M-2)]
				*vhx[i-1+j*M+(P-1)*N*M]+rescoefpp3boundF[i-1+(j-1)*(M-2)]*vhx[i+j*M+M*N]+rescoefpm3boundF[i-1+(j-1)*(M-2)]*vhx[i+j*M+(P-2)*M*N]
				+viscpp2[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i+1+(j+1)*M+(P-1)*N*M]-vhy[i+1+(j-1)*M+(P-1)*N*M])+viscpm2[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i-1+(j-1)*M+(P-1)*N*M]
				-vhy[i-1+(j+1)*M+(P-1)*N*M])+viscpp3[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhz[i+(j+1)*M+M*N]-vhz[i+(j-1)*M+M*N])+viscpm3[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]
				*(vhz[i+(j-1)*M+(P-2)*M*N]-vhz[i+(j+1)*M+(P-2)*M*N]))/(coefxboundF[i-1+(j-1)*(M-2)]);

		vhx[i+M*j]=vhx[i+M*j+M*N*(P-1)];
		}
	}
#pragma omp for
for(i=2; i<=M-3; i=i+2){             
	vhx[i+M*(N-1)+N*M*(P-1)]=(fhx[i+M*(N-1)+M*N*(P-1)]+2*rescoefpp1boundFR[i-1]*vhx[i+M+(P-1)*M*N]+
		2*rescoefpm1boundFR[i-1]*vhx[i+(N-2)*M+(P-1)*M*N]+rescoefpp2boundFR[i-1]*vhx[i+1+(N-1)*M+(P-1)*N*M]+rescoefpm2boundFR[i-1]
		*vhx[i-1+(N-1)*M+(P-1)*N*M]+rescoefpp3boundFR[i-1]*vhx[i+(N-1)*M+M*N]+rescoefpm3boundFR[i-1]*vhx[i+(N-1)*M+(P-2)*M*N]
		+viscpp2[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i+1+M+(P-1)*N*M]-vhy[i+1+(N-2)*M+(P-1)*N*M])+viscpm2[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i-1+(N-2)*M+(P-1)*N*M]
		-vhy[i-1+M+(P-1)*N*M])+viscpp3[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhz[i+M+M*N]-vhz[i+(N-2)*M+M*N])+viscpm3[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]
		*(vhz[i+(N-2)*M+(P-2)*M*N]-vhz[i+M+(P-2)*M*N]))/(coefxboundFR[i-1]);
	
	vhx[i]=vhx[i+M*(N-1)+N*M*(P-1)];
	vhx[i+M*(N-1)]=vhx[i+M*(N-1)+N*M*(P-1)];
	vhx[i+M*N*(P-1)]=vhx[i+M*(N-1)+N*M*(P-1)];

	}

#pragma omp for private(i,j)
for(k=1; k<=P-2; k=k+2){
	for(j=2; j<=N-3; j=j+2){
		for (i=1; i<=M-2; i=i+2){
			vhy[i+j*M+k*N*M]=(fhy[i+j*M+k*N*M]+rescoefpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+(j+1)*M+k*M*N]+
				rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+(j-1)*M+k*M*N]+2*rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+1+j*M+k*N*M]+2*rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhy[i-1+j*M+k*M*N]+rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+j*M+(k+1)*M*N]+rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+j*M+(k-1)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+1+(j+1)*M+k*N*M]-vhx[i-1+(j+1)*M+k*N*M])+viscpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i-1+(j-1)*M+k*N*M]
				-vhx[i+1+(j-1)*M+k*N*M])+viscpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+1+j*M+(k+1)*M*N]-vhz[i-1+j*M+(k+1)*M*N])+viscpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhz[i-1+j*M+(k-1)*M*N]-vhz[i+1+j*M+(k-1)*M*N]))/(coefy[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);
				}
			}
		}
#pragma omp for private(i)
for(k=1; k<=P-2; k=k+2){
	for(i=1; i<=M-2; i=i+2){
		vhy[i+(N-1)*M+k*N*M]=(fhy[i+(N-1)*M+k*N*M]+rescoefpp1boundR[i-1+(k-1)*(M-2)]*vhy[i+M+k*M*N]+
			rescoefpm1boundR[i-1+(k-1)*(M-2)]*vhy[i+(N-2)*M+k*M*N]+2*rescoefpp2boundR[i-1+(k-1)*(M-2)]*vhy[i+1+(N-1)*M+k*N*M]+2*rescoefpm2boundR[i-1+(k-1)*(M-2)]				
			*vhy[i-1+(N-1)*M+k*M*N]+rescoefpp3boundR[i-1+(k-1)*(M-2)]*vhy[i+(N-1)*M+(k+1)*M*N]+rescoefpm3boundR[i-1+(k-1)*(M-2)]*vhy[i+(N-1)*M+(k-1)*M*N]
			+viscpp1[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+1+M+k*N*M]-vhx[i-1+M+k*N*M])+viscpm1[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i-1+(N-2)*M+k*N*M]
			-vhx[i+1+(N-2)*M+k*N*M])+viscpp3[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+1+(N-1)*M+(k+1)*M*N]-vhz[i-1+(N-1)*M+(k+1)*M*N])+viscpm3[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]
			*(vhz[i-1+(N-1)*M+(k-1)*M*N]-vhz[i+1+(N-1)*M+(k-1)*M*N]))/(coefyboundR[i-1+(k-1)*(M-2)]);
		
		vhy[i+N*M*k]=vhy[i+M*(N-1)+N*M*k];
		}
	}
#pragma omp for private(i,j)
for(k=1; k<=P-2; k=k+2){
	for(j=1; j<=N-2; j=j+2){
		for (i=2; i<=M-3; i=i+2){
			vhy[i+j*M+k*N*M]=(fhy[i+j*M+k*N*M]+rescoefpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+(j+1)*M+k*M*N]+
				rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+(j-1)*M+k*M*N]+2*rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+1+j*M+k*N*M]+2*rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhy[i-1+j*M+k*M*N]+rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+j*M+(k+1)*M*N]+rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+j*M+(k-1)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+1+(j+1)*M+k*N*M]-vhx[i-1+(j+1)*M+k*N*M])+viscpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i-1+(j-1)*M+k*N*M]
				-vhx[i+1+(j-1)*M+k*N*M])+viscpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+1+j*M+(k+1)*M*N]-vhz[i-1+j*M+(k+1)*M*N])+viscpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhz[i-1+j*M+(k-1)*M*N]-vhz[i+1+j*M+(k-1)*M*N]))/(coefy[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);
				}
			}
		}

#pragma omp for private(i,j)
for(k=2; k<=P-3; k=k+2){
	for(j=1; j<=N-2; j=j+2){
		for (i=1; i<=M-2; i=i+2){
			vhy[i+j*M+k*N*M]=(fhy[i+j*M+k*N*M]+rescoefpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+(j+1)*M+k*M*N]+
				rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+(j-1)*M+k*M*N]+2*rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+1+j*M+k*N*M]+2*rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhy[i-1+j*M+k*M*N]+rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+j*M+(k+1)*M*N]+rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+j*M+(k-1)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+1+(j+1)*M+k*N*M]-vhx[i-1+(j+1)*M+k*N*M])+viscpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i-1+(j-1)*M+k*N*M]
				-vhx[i+1+(j-1)*M+k*N*M])+viscpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+1+j*M+(k+1)*M*N]-vhz[i-1+j*M+(k+1)*M*N])+viscpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhz[i-1+j*M+(k-1)*M*N]-vhz[i+1+j*M+(k-1)*M*N]))/(coefy[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);
				}
			}
		}
#pragma omp for private(i)
for(j=1; j<=N-2; j=j+2){
	for(i=1; i<=M-2; i=i+2){
		vhy[i+j*M+(P-1)*N*M]=(fhy[i+j*M+(P-1)*N*M]+rescoefpp1boundF[i-1+(j-1)*(M-2)]*vhy[i+(j+1)*M+(P-1)*M*N]+
				rescoefpm1boundF[i-1+(j-1)*(M-2)]*vhy[i+(j-1)*M+(P-1)*M*N]+2*rescoefpp2boundF[i-1+(j-1)*(M-2)]*vhy[i+1+j*M+(P-1)*N*M]+2*rescoefpm2boundF[i-1+(j-1)*(M-2)]
				*vhy[i-1+j*M+(P-1)*M*N]+rescoefpp3boundF[i-1+(j-1)*(M-2)]*vhy[i+j*M+M*N]+rescoefpm3boundF[i-1+(j-1)*(M-2)]*vhy[i+j*M+(P-2)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+1+(j+1)*M+(P-1)*N*M]-vhx[i-1+(j+1)*M+(P-1)*N*M])+viscpm1[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i-1+(j-1)*M+(P-1)*N*M]
				-vhx[i+1+(j-1)*M+(P-1)*N*M])+viscpp3[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhz[i+1+j*M+M*N]-vhz[i-1+j*M+M*N])+viscpm3[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]
				*(vhz[i-1+j*M+(P-2)*M*N]-vhz[i+1+j*M+(P-2)*M*N]))/(coefyboundF[i-1+(j-1)*(M-2)]);
				
		vhy[i+M*j]=vhy[i+M*j+M*N*(P-1)];
		}
	}


#pragma omp for private(i,j)
for(k=2; k<=P-3; k=k+2){
	for(j=2; j<=N-3; j=j+2){
		for (i=2; i<=M-3; i=i+2){
			vhy[i+j*M+k*N*M]=(fhy[i+j*M+k*N*M]+rescoefpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+(j+1)*M+k*M*N]+
				rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+(j-1)*M+k*M*N]+2*rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+1+j*M+k*N*M]+2*rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhy[i-1+j*M+k*M*N]+rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+j*M+(k+1)*M*N]+rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+j*M+(k-1)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+1+(j+1)*M+k*N*M]-vhx[i-1+(j+1)*M+k*N*M])+viscpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i-1+(j-1)*M+k*N*M]
				-vhx[i+1+(j-1)*M+k*N*M])+viscpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+1+j*M+(k+1)*M*N]-vhz[i-1+j*M+(k+1)*M*N])+viscpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhz[i-1+j*M+(k-1)*M*N]-vhz[i+1+j*M+(k-1)*M*N]))/(coefy[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);
				}
			}
		}


#pragma omp for private(i)
for(k=2; k<=P-3; k=k+2){
	for(i=2; i<=M-3; i=i+2){
		vhy[i+(N-1)*M+k*N*M]=(fhy[i+(N-1)*M+k*N*M]+rescoefpp1boundR[i-1+(k-1)*(M-2)]*vhy[i+M+k*M*N]+
			rescoefpm1boundR[i-1+(k-1)*(M-2)]*vhy[i+(N-2)*M+k*M*N]+2*rescoefpp2boundR[i-1+(k-1)*(M-2)]*vhy[i+1+(N-1)*M+k*N*M]+2*rescoefpm2boundR[i-1+(k-1)*(M-2)]				
			*vhy[i-1+(N-1)*M+k*M*N]+rescoefpp3boundR[i-1+(k-1)*(M-2)]*vhy[i+(N-1)*M+(k+1)*M*N]+rescoefpm3boundR[i-1+(k-1)*(M-2)]*vhy[i+(N-1)*M+(k-1)*M*N]
			+viscpp1[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+1+M+k*N*M]-vhx[i-1+M+k*N*M])+viscpm1[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i-1+(N-2)*M+k*N*M]
			-vhx[i+1+(N-2)*M+k*N*M])+viscpp3[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+1+(N-1)*M+(k+1)*M*N]-vhz[i-1+(N-1)*M+(k+1)*M*N])+viscpm3[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]
			*(vhz[i-1+(N-1)*M+(k-1)*M*N]-vhz[i+1+(N-1)*M+(k-1)*M*N]))/(coefyboundR[i-1+(k-1)*(M-2)]);
		
		vhy[i+N*M*k]=vhy[i+M*(N-1)+N*M*k];
		}
	}
#pragma omp for private(i)
for(j=2; j<=N-3; j=j+2){
	for(i=2; i<=M-3; i=i+2){
		vhy[i+j*M+(P-1)*N*M]=(fhy[i+j*M+(P-1)*N*M]+rescoefpp1boundF[i-1+(j-1)*(M-2)]*vhy[i+(j+1)*M+(P-1)*M*N]+
				rescoefpm1boundF[i-1+(j-1)*(M-2)]*vhy[i+(j-1)*M+(P-1)*M*N]+2*rescoefpp2boundF[i-1+(j-1)*(M-2)]*vhy[i+1+j*M+(P-1)*N*M]+2*rescoefpm2boundF[i-1+(j-1)*(M-2)]
				*vhy[i-1+j*M+(P-1)*M*N]+rescoefpp3boundF[i-1+(j-1)*(M-2)]*vhy[i+j*M+M*N]+rescoefpm3boundF[i-1+(j-1)*(M-2)]*vhy[i+j*M+(P-2)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+1+(j+1)*M+(P-1)*N*M]-vhx[i-1+(j+1)*M+(P-1)*N*M])+viscpm1[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i-1+(j-1)*M+(P-1)*N*M]
				-vhx[i+1+(j-1)*M+(P-1)*N*M])+viscpp3[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhz[i+1+j*M+M*N]-vhz[i-1+j*M+M*N])+viscpm3[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]
				*(vhz[i-1+j*M+(P-2)*M*N]-vhz[i+1+j*M+(P-2)*M*N]))/(coefyboundF[i-1+(j-1)*(M-2)]);
				
		vhy[i+M*j]=vhy[i+M*j+M*N*(P-1)];
		}
	}
#pragma omp for
for(i=2; i<=M-3; i=i+2){
	vhy[i+(N-1)*M+(P-1)*N*M]=(fhy[i+(N-1)*M+(P-1)*N*M]+rescoefpp1boundFR[i-1]*vhy[i+M+(P-1)*M*N]+
		rescoefpm1boundFR[i-1]*vhy[i+(N-2)*M+(P-1)*M*N]+2*rescoefpp2boundFR[i-1]*vhy[i+1+(N-1)*M+(P-1)*N*M]+2*rescoefpm2boundFR[i-1]				
		*vhy[i-1+(N-1)*M+(P-1)*M*N]+rescoefpp3boundFR[i-1]*vhy[i+(N-1)*M+M*N]+rescoefpm3boundFR[i-1]*vhy[i+(N-1)*M+(P-2)*M*N]
		+viscpp1[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+1+M+(P-1)*N*M]-vhx[i-1+M+(P-1)*N*M])+viscpm1[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i-1+(N-2)*M+(P-1)*N*M]
		-vhx[i+1+(N-2)*M+(P-1)*N*M])+viscpp3[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhz[i+1+(N-1)*M+M*N]-vhz[i-1+(N-1)*M+M*N])+viscpm3[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]
		*(vhz[i-1+(N-1)*M+(P-2)*M*N]-vhz[i+1+(N-1)*M+(P-2)*M*N]))/(coefyboundFR[i-1]);

	vhy[i]=vhy[i+M*(N-1)+N*M*(P-1)];
	vhy[i+M*(N-1)]=vhy[i+M*(N-1)+N*M*(P-1)];
	vhy[i+M*N*(P-1)]=vhy[i+M*(N-1)+N*M*(P-1)];
	}
#pragma omp for private(i,j)
for(k=1; k<=P-2; k=k+2){
	for(j=2; j<=N-3; j=j+2){
		for(i=1; i<=M-2; i=i+2){
			vhz[i+j*M+k*N*M]=(fhz[i+j*M+k*N*M]+rescoefpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+(j+1)*M+k*M*N]+
				rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+(j-1)*M+k*M*N]+rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+1+j*M+k*N*M]+rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhz[i-1+j*M+k*M*N]+2*rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+j*M+(k+1)*M*N]+2*rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+j*M+(k-1)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(j+1)*M+(k+1)*N*M]-vhx[i+(j+1)*M+(k-1)*N*M])+viscpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(j-1)*M+(k-1)*N*M]
				-vhx[i+(j-1)*M+(k+1)*N*M])+viscpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+j*M+(k+1)*M*N]-vhy[i+1+j*M+(k-1)*M*N])+viscpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhy[i-1+j*M+(k-1)*M*N]-vhy[i-1+j*M+(k+1)*M*N]))/(coefz[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);

			}
		}
	}
#pragma omp for private(i)
for(k=1; k<=P-2; k=k+2){
	for(i=1; i<=M-2; i=i+2){
		vhz[i+(N-1)*M+k*N*M]=(fhz[i+(N-1)*M+k*N*M]+rescoefpp1boundR[i-1+(k-1)*(M-2)]*vhz[i+M+k*M*N]+
			rescoefpm1boundR[i-1+(k-1)*(M-2)]*vhz[i+(N-2)*M+k*M*N]+rescoefpp2boundR[i-1+(k-1)*(M-2)]*vhz[i+1+(N-1)*M+k*N*M]+rescoefpm2boundR[i-1+(k-1)*(M-2)]
			*vhz[i-1+(N-1)*M+k*M*N]+2*rescoefpp3boundR[i-1+(k-1)*(M-2)]*vhz[i+(N-1)*M+(k+1)*M*N]+2*rescoefpm3boundR[i-1+(k-1)*(M-2)]*vhz[i+(N-1)*M+(k-1)*M*N]
			+viscpp1[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+M+(k+1)*N*M]-vhx[i+M+(k-1)*N*M])+viscpm1[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(N-2)*M+(k-1)*N*M]
			-vhx[i+(N-2)*M+(k+1)*N*M])+viscpp2[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+(N-1)*M+(k+1)*M*N]-vhy[i+1+(N-1)*M+(k-1)*M*N])+viscpm2[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]
			*(vhy[i-1+(N-1)*M+(k-1)*M*N]-vhy[i-1+(N-1)*M+(k+1)*M*N]))/(coefzboundR[i-1+(k-1)*(M-2)]);

		vhz[i+N*M*k]=vhz[i+M*(N-1)+N*M*k];
		}
	}
#pragma omp for private(i,j)
for(k=1; k<=P-2; k=k+2){
	for(j=1; j<=N-2; j=j+2){
		for(i=2; i<=M-3; i=i+2){
			vhz[i+j*M+k*N*M]=(fhz[i+j*M+k*N*M]+rescoefpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+(j+1)*M+k*M*N]+
				rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+(j-1)*M+k*M*N]+rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+1+j*M+k*N*M]+rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhz[i-1+j*M+k*M*N]+2*rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+j*M+(k+1)*M*N]+2*rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+j*M+(k-1)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(j+1)*M+(k+1)*N*M]-vhx[i+(j+1)*M+(k-1)*N*M])+viscpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(j-1)*M+(k-1)*N*M]
				-vhx[i+(j-1)*M+(k+1)*N*M])+viscpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+j*M+(k+1)*M*N]-vhy[i+1+j*M+(k-1)*M*N])+viscpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhy[i-1+j*M+(k-1)*M*N]-vhy[i-1+j*M+(k+1)*M*N]))/(coefz[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);

			}
		}
	}
#pragma omp for private(i,j)
for(k=2; k<=P-3; k=k+2){
	for(j=1; j<=N-2; j=j+2){
		for(i=1; i<=M-2; i=i+2){
			vhz[i+j*M+k*N*M]=(fhz[i+j*M+k*N*M]+rescoefpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+(j+1)*M+k*M*N]+
				rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+(j-1)*M+k*M*N]+rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+1+j*M+k*N*M]+rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhz[i-1+j*M+k*M*N]+2*rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+j*M+(k+1)*M*N]+2*rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+j*M+(k-1)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(j+1)*M+(k+1)*N*M]-vhx[i+(j+1)*M+(k-1)*N*M])+viscpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(j-1)*M+(k-1)*N*M]
				-vhx[i+(j-1)*M+(k+1)*N*M])+viscpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+j*M+(k+1)*M*N]-vhy[i+1+j*M+(k-1)*M*N])+viscpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhy[i-1+j*M+(k-1)*M*N]-vhy[i-1+j*M+(k+1)*M*N]))/(coefz[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);

			}
		}
	}


#pragma omp for private(i)
for(j=1; j<=N-2; j=j+2){
	for(i=1; i<=M-2; i=i+2){
		vhz[i+j*M+(P-1)*N*M]=(fhz[i+j*M+(P-1)*N*M]+rescoefpp1boundF[i-1+(j-1)*(M-2)]*vhz[i+(j+1)*M+(P-1)*M*N]+
				rescoefpm1boundF[i-1+(j-1)*(M-2)]*vhz[i+(j-1)*M+(P-1)*M*N]+rescoefpp2boundF[i-1+(j-1)*(M-2)]*vhz[i+1+j*M+(P-1)*N*M]+rescoefpm2boundF[i-1+(j-1)*(M-2)]
				*vhz[i-1+j*M+(P-1)*M*N]+2*rescoefpp3boundF[i-1+(j-1)*(M-2)]*vhz[i+j*M+M*N]+2*rescoefpm3boundF[i-1+(j-1)*(M-2)]*vhz[i+j*M+M*N*(P-2)]
				+viscpp1[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+(j+1)*M+M*N]-vhx[i+(j+1)*M+(P-2)*N*M])+viscpm1[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+(j-1)*M+(P-2)*N*M]
				-vhx[i+(j-1)*M+M*N])+viscpp2[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i+1+j*M+M*N]-vhy[i+1+j*M+(P-2)*M*N])+viscpm2[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]
				*(vhy[i-1+j*M+(P-2)*M*N]-vhy[i-1+j*M+M*N]))/(coefzboundF[i-1+(j-1)*(M-2)]);

		vhz[i+M*j]=vhz[i+M*j+M*N*(P-1)];
		}
	}
#pragma omp for private(i,j)
for(k=2; k<=P-3; k=k+2){
	for(j=2; j<=N-3; j=j+2){
		for(i=2; i<=M-3; i=i+2){
			vhz[i+j*M+k*N*M]=(fhz[i+j*M+k*N*M]+rescoefpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+(j+1)*M+k*M*N]+
				rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+(j-1)*M+k*M*N]+rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+1+j*M+k*N*M]+rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhz[i-1+j*M+k*M*N]+2*rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+j*M+(k+1)*M*N]+2*rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+j*M+(k-1)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(j+1)*M+(k+1)*N*M]-vhx[i+(j+1)*M+(k-1)*N*M])+viscpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(j-1)*M+(k-1)*N*M]
				-vhx[i+(j-1)*M+(k+1)*N*M])+viscpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+j*M+(k+1)*M*N]-vhy[i+1+j*M+(k-1)*M*N])+viscpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhy[i-1+j*M+(k-1)*M*N]-vhy[i-1+j*M+(k+1)*M*N]))/(coefz[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]);

			}
		}
	}
#pragma omp for private(i)
for(k=2; k<=P-3; k=k+2){
	for(i=2; i<=M-3; i=i+2){
		vhz[i+(N-1)*M+k*N*M]=(fhz[i+(N-1)*M+k*N*M]+rescoefpp1boundR[i-1+(k-1)*(M-2)]*vhz[i+M+k*M*N]+
			rescoefpm1boundR[i-1+(k-1)*(M-2)]*vhz[i+(N-2)*M+k*M*N]+rescoefpp2boundR[i-1+(k-1)*(M-2)]*vhz[i+1+(N-1)*M+k*N*M]+rescoefpm2boundR[i-1+(k-1)*(M-2)]
			*vhz[i-1+(N-1)*M+k*M*N]+2*rescoefpp3boundR[i-1+(k-1)*(M-2)]*vhz[i+(N-1)*M+(k+1)*M*N]+2*rescoefpm3boundR[i-1+(k-1)*(M-2)]*vhz[i+(N-1)*M+(k-1)*M*N]
			+viscpp1[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+M+(k+1)*N*M]-vhx[i+M+(k-1)*N*M])+viscpm1[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(N-2)*M+(k-1)*N*M]
			-vhx[i+(N-2)*M+(k+1)*N*M])+viscpp2[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+(N-1)*M+(k+1)*M*N]-vhy[i+1+(N-1)*M+(k-1)*M*N])+viscpm2[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]
			*(vhy[i-1+(N-1)*M+(k-1)*M*N]-vhy[i-1+(N-1)*M+(k+1)*M*N]))/(coefzboundR[i-1+(k-1)*(M-2)]);

		vhz[i+N*M*k]=vhz[i+M*(N-1)+N*M*k];
		}
	}
#pragma omp for private(i)
for(j=2; j<=N-3; j=j+2){
	for(i=2; i<=M-3; i=i+2){
		vhz[i+j*M+(P-1)*N*M]=(fhz[i+j*M+(P-1)*N*M]+rescoefpp1boundF[i-1+(j-1)*(M-2)]*vhz[i+(j+1)*M+(P-1)*M*N]+
				rescoefpm1boundF[i-1+(j-1)*(M-2)]*vhz[i+(j-1)*M+(P-1)*M*N]+rescoefpp2boundF[i-1+(j-1)*(M-2)]*vhz[i+1+j*M+(P-1)*N*M]+rescoefpm2boundF[i-1+(j-1)*(M-2)]
				*vhz[i-1+j*M+(P-1)*M*N]+2*rescoefpp3boundF[i-1+(j-1)*(M-2)]*vhz[i+j*M+M*N]+2*rescoefpm3boundF[i-1+(j-1)*(M-2)]*vhz[i+j*M+M*N*(P-2)]
				+viscpp1[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+(j+1)*M+M*N]-vhx[i+(j+1)*M+(P-2)*N*M])+viscpm1[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+(j-1)*M+(P-2)*N*M]
				-vhx[i+(j-1)*M+M*N])+viscpp2[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i+1+j*M+M*N]-vhy[i+1+j*M+(P-2)*M*N])+viscpm2[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]
				*(vhy[i-1+j*M+(P-2)*M*N]-vhy[i-1+j*M+M*N]))/(coefzboundF[i-1+(j-1)*(M-2)]);

		vhz[i+M*j]=vhz[i+M*j+M*N*(P-1)];
		}
	}
#pragma omp for
for(i=2; i<=M-3; i=i+2){
	vhz[i+(N-1)*M+(P-1)*N*M]=(fhz[i+(N-1)*M+(P-1)*N*M]+rescoefpp1boundFR[i-1]*vhz[i+M+(P-1)*M*N]+
		rescoefpm1boundFR[i-1]*vhz[i+(N-2)*M+(P-1)*M*N]+rescoefpp2boundFR[i-1]*vhz[i+1+(N-1)*M+(P-1)*N*M]+rescoefpm2boundFR[i-1]
		*vhz[i-1+(N-1)*M+(P-1)*M*N]+2*rescoefpp3boundFR[i-1]*vhz[i+(N-1)*M+M*N]+2*rescoefpm3boundFR[i-1]*vhz[i+(N-1)*M+(P-2)*M*N]
		+viscpp1[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+M+M*N]-vhx[i+M+(P-2)*N*M])+viscpm1[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+(N-2)*M+(P-2)*N*M]
		-vhx[i+(N-2)*M+M*N])+viscpp2[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i+1+(N-1)*M+M*N]-vhy[i+1+(N-1)*M+(P-2)*M*N])+viscpm2[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]
		*(vhy[i-1+(N-1)*M+(P-2)*M*N]-vhy[i-1+(N-1)*M+M*N]))/(coefzboundFR[i-1]);

	vhz[i]=vhz[i+M*(N-1)+N*M*(P-1)];
	vhz[i+M*(N-1)]=vhz[i+M*(N-1)+N*M*(P-1)];
	vhz[i+M*N*(P-1)]=vhz[i+M*(N-1)+N*M*(P-1)];

	}

}

}

plhs[0]=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS, mxREAL);
plhs[1]=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS, mxREAL);
plhs[2]=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS, mxREAL);

memcpy(mxGetPr(plhs[0]),vhx,sizeof(double)*M*N*P);
memcpy(mxGetPr(plhs[1]),vhy,sizeof(double)*M*N*P);
memcpy(mxGetPr(plhs[2]),vhz,sizeof(double)*M*N*P);

return;
}
