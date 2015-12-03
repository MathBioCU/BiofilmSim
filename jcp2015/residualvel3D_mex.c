#include <math.h>
#include <string.h> /*for memcpy function*/
#include "mex.h"
#include "matrix.h"
#include <omp.h>


 #define min(X,Y) ((X) < (Y) ? (X) : (Y))
/*#define vhx(i,j,k) vhx[(i)+M*(j)+N*M*(k)]
#define vhy(i,j,k) vhy[(i)+M*(j)+N*M*(k)]
#define vhz(i,j,k) vhz[(i)+M*(j)+N*M*(k)]
#define fhx(i,j,k) fhx[(i)+M*(j)+N*M*(k)]
#define fhy(i,j,k) fhy[(i)+M*(j)+N*M*(k)]
#define fhz(i,j,k) fhz[(i)+M*(j)+N*M*(k)]
#define rx(i,j,k) rx[(i)+M*(j)+N*M*(k)]
#define ry(i,j,k) ry[(i)+M*(j)+N*M*(k)]
#define rz(i,j,k) rz[(i)+M*(j)+N*M*(k)]*/
/*#define (A)(i,j,k,M,N) (A)[(i)+(M)*(j)+(N)*(M)*(k)]*/
/*#define coefx(i,j,k,M,N) coefx[(i)+(M)*((j))+(M)*(N)*((k))]
#define coefy(i,j,k,M,N) coefy[(i)+(M)*((j))+(M)*(N)*((k))]
#define coefz(i,j,k,M,N) coefz[(i)+(M)*((j))+(M)*(N)*((k))]
*/
/*need to be careful not to define macros that overwrite coef vs. coefbound etc.*/
/*#define rescoefpp1(i,j,k,M,N) rescoefpp1[(i)+(M)*((j))+(M)*(N)*((k))]
#define rescoefpm1(i,j,k,M,N) rescoefpm1[(i)+(M)*((j))+(M)*(N)*((k))]
#define rescoefpp2(i,j,k,M,N) rescoefpp2[(i)+(M)*((j))+(M)*(N)*((k))]
#define rescoefpm2(i,j,k,M,N) rescoefpm2[(i)+(M)*((j))+(M)*(N)*((k))]
#define rescoefpp3(i,j,k,M,N) rescoefpp3[(i)+(M)*((j))+(M)*(N)*((k))]
#define rescoefpm3(i,j,k,M,N) rescoefpm3[(i)+(M)*((j))+(M)*(N)*((k))]
#define viscpp1(i,j,k,M,N)  viscpp1[(i)+(M)*(j)+(M)*(N)*(k)]
#define viscpm1(i,j,k,M,N)  viscpm1[(i)+(M)*(j)+(M)*(N)*(k)]
N(X,Y) ((X) < (Y) ? : (X) : (Y))#define viscpp2(i,j,k,M,N)  viscpp2[(i)+(M)*(j)+(M)*(N)*(k)]
N(X,Y) ((X) < (Y) ? : (X) : (Y))N(X,Y) ((X) < (Y) ? (X) : (Y))#define viscpm2(i,j,k,M,N)  viscpm2[(i)+(M)*(j)+(M)*(N)*(k)]
#define viscpp3(i,j,k,M,N)  viscpp3[(i)+(M)*(j)+(M)*(N)*(k)]
# viscpm3(i,j,k,M,N)  viscpm3[(i)+(M)*(j)+(M)*(N)*(k)]
*/

/*
struct coef
{
double *coefx, *coefy, *coefz, *rescoefpm1, *rescoefpm2, *rescoefpm3, *rescoefpp1;
double *rescoefpp2, *rescoefpp3, *viscpp1, *viscpp2, *viscpp3, *viscpm1, *viscpm2; 
double *viscpm3, *rescoefpm1boundF, *rescoefpp1boundF, *rescoefpp2boundF;
double *rescoefpm2boundF, *rescoefpp3boundF, *rescoefpm3boundF, *rescoefpp1boundR;
double *rescoefpm1boundR, *rescoefpm2boundR, *rescoefpp2boundR, *rescoefpm3boundR;
double *rescoefpp3boundR, *coefxboundR, *coefyboundR, *coefzboundR, *coefxboundF;
double *rescoefpp1boundFR, *rescoefpm1boundFR, *rescoefpp2boundFR, *rescoefpm2boundFR;
double *rescoefpp3boundFR, *rescoefpm3boundFR;
double *coefyboundF, *coefzboundF, *coefxboundFR, *coefyboundFR, *coefzboundFR; 
};*/

void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
int i,j,k;
double *vhx, *vhy,*vhz,*fhx,*fhy,*fhz, *rx, *ry, *rz;
double *coefx,*coefy,*coefz, *rescoefpm1, *rescoefpm2, *rescoefpm3, *rescoefpp1, *rescoefpp2, *rescoefpp3, *viscpp1, *viscpp2, *viscpp3, *viscpm1, *viscpm2, *viscpm3;
double *rescoefpm1boundF, *rescoefpp1boundF,*rescoefpp2boundF, *rescoefpm2boundF, *rescoefpp3boundF, *rescoefpm3boundF, *rescoefpp1boundR, *rescoefpm1boundR;
double *rescoefpm2boundR, *rescoefpp2boundR, *rescoefpm3boundR, *rescoefpp3boundR, *coefxboundR, *coefyboundR, *coefzboundR, *coefxboundF, *coefyboundF, *coefzboundF;
double *coefxboundFR, *coefyboundFR, *coefzboundFR, *rescoefpp1boundFR, *rescoefpp2boundFR, *rescoefpp3boundFR, *rescoefpm1boundFR, *rescoefpm2boundFR, *rescoefpm3boundFR;


/* Input Arguments */
vhx = mxGetPr(prhs[0]);
vhy = mxGetPr(prhs[1]);
vhz = mxGetPr(prhs[2]);
fhx = mxGetPr(prhs[3]);
fhy = mxGetPr(prhs[4]);
fhz = mxGetPr(prhs[5]);
rx = mxGetPr(prhs[7]);
ry = mxGetPr(prhs[8]);
rz = mxGetPr(prhs[9]);

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

int r,s,t, st1, st2, st3, en1, en2, en3, bs;
bs=25;
if(omp_get_num_procs()>=4){
	omp_set_num_threads(4);}
#pragma omp parallel
{
#pragma omp for private(i,j)
for(r=1; r<=P-2; r=r+bs){
        st1=r; 
	en1=min(r+bs-1, P-2);
	for(s=1; s<=N-2; s+bs){
	     st2=s; en2=min(s+bs-1, N-2);
		for (t=1; t<=M-2; t+bs){
		     st3=t; en3=min(t+bs-1, M-2);
			for(i=st1; i<en1; i++){
				for(j=st2; j<en2; j++){
					for(k=st3; k<en3; k++){ 
		


			rx[i+M*j+N*M*k]=fhx[i+M*j+M*N*k]-coefx[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+M*j+M*N*k]+2*rescoefpp1[i-1+(M-2)*(j-1)+(k-1)*(M-2)*(N-2)]*vhx[i+(j+1)*M+k*M*N]+
				2*rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+(j-1)*M+k*M*N]+rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+1+j*M+k*N*M]+rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhx[i-1+j*M+k*N*M]+rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+j*M+(k+1)*M*N]+rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhx[i+j*M+(k-1)*M*N]
				+viscpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+(j+1)*M+k*N*M]-vhy[i+1+(j-1)*M+k*N*M])+viscpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i-1+(j-1)*M+k*N*M]
				-vhy[i-1+(j+1)*M+k*N*M])+viscpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+(j+1)*M+(k+1)*M*N]-vhz[i+(j-1)*M+(k+1)*M*N])+viscpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhz[i+(j-1)*M+(k-1)*M*N]-vhz[i+(j+1)*M+(k-1)*M*N]);
							}
						}
					}
				}
			}
		}


#pragma omp for private(i,j)
for(r=1; r<=P-2; r=r+bs){
        st1=r; en1=min(r+bs-1, P-2);
	for(s=1; s<=N-2; s+bs){
	     st2=s; en2=min(s+bs-1, N-2);
		for (t=1; t<=M-2; t+bs){
		     st3=t; en3=min(t+bs-1, M-2);
			for(i=st1; i<en1; i++){
				for(j=st2; j<en2; j++){
					for(k=st3; k<en3; k++){ 

			ry[i+j*M+k*N*M]=fhy[i+j*M+k*N*M]-coefy[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+j*M+k*N*M]+rescoefpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+(j+1)*M+k*M*N]+
				rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+(j-1)*M+k*M*N]+2*rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+1+j*M+k*N*M]+2*rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhy[i-1+j*M+k*M*N]+rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+j*M+(k+1)*M*N]+rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhy[i+j*M+(k-1)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+1+(j+1)*M+k*N*M]-vhx[i-1+(j+1)*M+k*N*M])+viscpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i-1+(j-1)*M+k*N*M]
				-vhx[i+1+(j-1)*M+k*N*M])+viscpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+1+j*M+(k+1)*M*N]-vhz[i-1+j*M+(k+1)*M*N])+viscpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhz[i-1+j*M+(k-1)*M*N]-vhz[i+1+j*M+(k-1)*M*N]);
							}
						}
					}
				}
			}
		}


#pragma omp for private(i,j)
for(r=1; r<=P-2; r=r+bs){
        st1=r; en1=min(r+bs-1, P-2);
	for(s=1; s<=N-2; s+bs){
	     st2=s; en2=min(s+bs-1, N-2);
		for (t=1; t<=M-2; t+bs){
		     st3=t; en3=min(t+bs-1, M-2);
			for(i=st1; i<en1; i++){
				for(j=st2; j<en2; j++){
					for(k=st3; k<en3; k++){ 

			rz[i+j*M+k*N*M]=fhz[i+j*M+k*N*M]-coefz[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+j*M+k*N*M]+rescoefpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+(j+1)*M+k*M*N]+
				rescoefpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+(j-1)*M+k*M*N]+rescoefpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+1+j*M+k*N*M]+rescoefpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]
				*vhz[i-1+j*M+k*M*N]+2*rescoefpp3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+j*M+(k+1)*M*N]+2*rescoefpm3[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-2)]*vhz[i+j*M+(k-1)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(j+1)*M+(k+1)*N*M]-vhx[i+(j+1)*M+(k-1)*N*M])+viscpm1[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(j-1)*M+(k-1)*N*M]
				-vhx[i+(j-1)*M+(k+1)*N*M])+viscpp2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+j*M+(k+1)*M*N]-vhy[i+1+j*M+(k-1)*M*N])+viscpm2[i-1+(j-1)*(M-2)+(k-1)*(M-2)*(N-1)]
				*(vhy[i-1+j*M+(k-1)*M*N]-vhy[i-1+j*M+(k+1)*M*N]);
						}
					}
				}
					
			}
		}
	}

}

for(j=1; j<=N-2; j++){
	for(i=1; i<=M-2; i++){
		

		rx[i+M*j+N*M*(P-1)]=fhx[i+M*j+M*N*(P-1)]-coefxboundF[i-1+(j-1)*(M-2)]*vhx[i+M*j+M*N*(P-1)]+2*rescoefpp1boundF[i-1+(M-2)*(j-1)]*vhx[i+(j+1)*M+(P-1)*M*N]+
				2*rescoefpm1boundF[i-1+(j-1)*(M-2)]*vhx[i+(j-1)*M+(P-1)*M*N]+rescoefpp2boundF[i-1+(j-1)*(M-2)]*vhx[i+1+j*M+(P-1)*N*M]+rescoefpm2boundF[i-1+(j-1)*(M-2)]
				*vhx[i-1+j*M+(P-1)*N*M]+rescoefpp3boundF[i-1+(j-1)*(M-2)]*vhx[i+j*M+M*N]+rescoefpm3boundF[i-1+(j-1)*(M-2)]*vhx[i+j*M+(P-2)*M*N]
				+viscpp2[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i+1+(j+1)*M+(P-1)*N*M]-vhy[i+1+(j-1)*M+(P-1)*N*M])+viscpm2[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i-1+(j-1)*M+(P-1)*N*M]
				-vhy[i-1+(j+1)*M+(P-1)*N*M])+viscpp3[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhz[i+(j+1)*M+M*N]-vhz[i+(j-1)*M+M*N])+viscpm3[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]
				*(vhz[i+(j-1)*M+(P-2)*M*N]-vhz[i+(j+1)*M+(P-2)*M*N]);

		rx[i+M*j]=rx[i+M*j+M*N*(P-1)];

		ry[i+j*M+(P-1)*N*M]=fhy[i+j*M+(P-1)*N*M]-coefyboundF[i-1+(j-1)*(M-2)]*vhy[i+j*M+(P-1)*N*M]+rescoefpp1boundF[i-1+(j-1)*(M-2)]*vhy[i+(j+1)*M+(P-1)*M*N]+
				rescoefpm1boundF[i-1+(j-1)*(M-2)]*vhy[i+(j-1)*M+(P-1)*M*N]+2*rescoefpp2boundF[i-1+(j-1)*(M-2)]*vhy[i+1+j*M+(P-1)*N*M]+2*rescoefpm2boundF[i-1+(j-1)*(M-2)]
				*vhy[i-1+j*M+(P-1)*M*N]+rescoefpp3boundF[i-1+(j-1)*(M-2)]*vhy[i+j*M+M*N]+rescoefpm3boundF[i-1+(j-1)*(M-2)]*vhy[i+j*M+(P-2)*M*N]
				+viscpp1[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+1+(j+1)*M+(P-1)*N*M]-vhx[i-1+(j+1)*M+(P-1)*N*M])+viscpm1[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i-1+(j-1)*M+(P-1)*N*M]
				-vhx[i+1+(j-1)*M+(P-1)*N*M])+viscpp3[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhz[i+1+j*M+M*N]-vhz[i-1+j*M+M*N])+viscpm3[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]
				*(vhz[i-1+j*M+(P-2)*M*N]-vhz[i+1+j*M+(P-2)*M*N]);
				
		
		ry[i+M*j]=ry[i+M*j+M*N*(P-1)];

		rz[i+j*M+(P-1)*N*M]=fhz[i+j*M+(P-1)*N*M]-coefzboundF[i-1+(j-1)*(M-2)]*vhz[i+j*M+(P-1)*N*M]+rescoefpp1boundF[i-1+(j-1)*(M-2)]*vhz[i+(j+1)*M+(P-1)*M*N]+
				rescoefpm1boundF[i-1+(j-1)*(M-2)]*vhz[i+(j-1)*M+(P-1)*M*N]+rescoefpp2boundF[i-1+(j-1)*(M-2)]*vhz[i+1+j*M+(P-1)*N*M]+rescoefpm2boundF[i-1+(j-1)*(M-2)]
				*vhz[i-1+j*M+(P-1)*M*N]+2*rescoefpp3boundF[i-1+(j-1)*(M-2)]*vhz[i+j*M+M*N]+2*rescoefpm3boundF[i-1+(j-1)*(M-2)]*vhz[i+j*M+M*N*(P-2)]
				+viscpp1[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+(j+1)*M+M*N]-vhx[i+(j+1)*M+(P-2)*N*M])+viscpm1[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+(j-1)*M+(P-2)*N*M]
				-vhx[i+(j-1)*M+M*N])+viscpp2[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i+1+j*M+M*N]-vhy[i+1+j*M+(P-2)*M*N])+viscpm2[i-1+(j-1)*(M-2)+(P-2)*(M-2)*(N-1)]
				*(vhy[i-1+j*M+(P-2)*M*N]-vhy[i-1+j*M+M*N]);

		rz[i+M*j]=rz[i+M*j+M*N*(P-1)];

		}
	}


for(k=1; k<=P-2; k++){
	for(i=1; i<=M-2; i++){

		rx[i+M*(N-1)+N*M*k]=fhx[i+M*(N-1)+M*N*k]-coefxboundR[i-1+(k-1)*(M-2)]*vhx[i+M*(N-1)+M*N*k]+2*rescoefpp1boundR[i-1+(k-1)*(M-2)]*vhx[i+M+k*M*N]+
			2*rescoefpm1boundR[i-1+(k-1)*(M-2)]*vhx[i+(N-2)*M+k*M*N]+rescoefpp2boundR[i-1+(k-1)*(M-2)]*vhx[i+1+(N-1)*M+k*N*M]+rescoefpm2boundR[i-1+(k-1)*(M-2)]
			*vhx[i-1+(N-1)*M+k*N*M]+rescoefpp3boundR[i-1+(k-1)*(M-2)]*vhx[i+(N-1)*M+(k+1)*M*N]+rescoefpm3boundR[i-1+(k-1)*(M-2)]*vhx[i+(N-1)*M+(k-1)*M*N]
			+viscpp2[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+M+k*N*M]-vhy[i+1+(N-2)*M+k*N*M])+viscpm2[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i-1+(N-2)*M+k*N*M]
			-vhy[i-1+M+k*N*M])+viscpp3[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+M+(k+1)*M*N]-vhz[i+(N-2)*M+(k+1)*M*N])+viscpm3[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]
			*(vhz[i+(N-2)*M+(k-1)*M*N]-vhz[i+M+(k-1)*M*N]);
		rx[i+N*M*k]=rx[i+M*(N-1)+N*M*k];


		ry[i+(N-1)*M+k*N*M]=fhy[i+(N-1)*M+k*N*M]-coefyboundR[i-1+(k-1)*(M-2)]*vhy[i+(N-1)*M+k*N*M]+rescoefpp1boundR[i-1+(k-1)*(M-2)]*vhy[i+M+k*M*N]+
			rescoefpm1boundR[i-1+(k-1)*(M-2)]*vhy[i+(N-2)*M+k*M*N]+2*rescoefpp2boundR[i-1+(k-1)*(M-2)]*vhy[i+1+(N-1)*M+k*N*M]+2*rescoefpm2boundR[i-1+(k-1)*(M-2)]				
			*vhy[i-1+(N-1)*M+k*M*N]+rescoefpp3boundR[i-1+(k-1)*(M-2)]*vhy[i+(N-1)*M+(k+1)*M*N]+rescoefpm3boundR[i-1+(k-1)*(M-2)]*vhy[i+(N-1)*M+(k-1)*M*N]
			+viscpp1[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+1+M+k*N*M]-vhx[i-1+M+k*N*M])+viscpm1[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i-1+(N-2)*M+k*N*M]
			-vhx[i+1+(N-2)*M+k*N*M])+viscpp3[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhz[i+1+(N-1)*M+(k+1)*M*N]-vhz[i-1+(N-1)*M+(k+1)*M*N])+viscpm3[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]
			*(vhz[i-1+(N-1)*M+(k-1)*M*N]-vhz[i+1+(N-1)*M+(k-1)*M*N]);
		
		ry[i+N*M*k]=ry[i+M*(N-1)+N*M*k];

		
		rz[i+(N-1)*M+k*N*M]=fhz[i+(N-1)*M+k*N*M]-coefzboundR[i-1+(k-1)*(M-2)]*vhz[i+(N-1)*M+k*N*M]+rescoefpp1boundR[i-1+(k-1)*(M-2)]*vhz[i+M+k*M*N]+
			rescoefpm1boundR[i-1+(k-1)*(M-2)]*vhz[i+(N-2)*M+k*M*N]+rescoefpp2boundR[i-1+(k-1)*(M-2)]*vhz[i+1+(N-1)*M+k*N*M]+rescoefpm2boundR[i-1+(k-1)*(M-2)]
			*vhz[i-1+(N-1)*M+k*M*N]+2*rescoefpp3boundR[i-1+(k-1)*(M-2)]*vhz[i+(N-1)*M+(k+1)*M*N]+2*rescoefpm3boundR[i-1+(k-1)*(M-2)]*vhz[i+(N-1)*M+(k-1)*M*N]
			+viscpp1[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+M+(k+1)*N*M]-vhx[i+M+(k-1)*N*M])+viscpm1[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhx[i+(N-2)*M+(k-1)*N*M]
			-vhx[i+(N-2)*M+(k+1)*N*M])+viscpp2[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]*(vhy[i+1+(N-1)*M+(k+1)*M*N]-vhy[i+1+(N-1)*M+(k-1)*M*N])+viscpm2[i-1+(N-2)*(M-2)+(k-1)*(M-2)*(N-1)]
			*(vhy[i-1+(N-1)*M+(k-1)*M*N]-vhy[i-1+(N-1)*M+(k+1)*M*N]);

		rz[i+N*M*k]=rz[i+M*(N-1)+N*M*k];

		}
	}




for(i=1; i<=M-2; i++){             
	rx[i+M*(N-1)+N*M*(P-1)]=fhx[i+M*(N-1)+M*N*(P-1)]-coefxboundFR[i-1]*vhx[i+M*(N-1)+M*N*(P-1)]+2*rescoefpp1boundFR[i-1]*vhx[i+M+(P-1)*M*N]+
		2*rescoefpm1boundFR[i-1]*vhx[i+(N-2)*M+(P-1)*M*N]+rescoefpp2boundFR[i-1]*vhx[i+1+(N-1)*M+(P-1)*N*M]+rescoefpm2boundFR[i-1]
		*vhx[i-1+(N-1)*M+(P-1)*N*M]+rescoefpp3boundFR[i-1]*vhx[i+(N-1)*M+M*N]+rescoefpm3boundFR[i-1]*vhx[i+(N-1)*M+(P-2)*M*N]
		+viscpp2[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i+1+M+(P-1)*N*M]-vhy[i+1+(N-2)*M+(P-1)*N*M])+viscpm2[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i-1+(N-2)*M+(P-1)*N*M]
		-vhy[i-1+M+(P-1)*N*M])+viscpp3[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhz[i+M+M*N]-vhz[i+(N-2)*M+M*N])+viscpm3[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]
		*(vhz[i+(N-2)*M+(P-2)*M*N]-vhz[i+M+(P-2)*M*N]);
	
	rx[i]=rx[i+M*(N-1)+N*M*(P-1)];
	rx[i+M*(N-1)]=rx[i+M*(N-1)+N*M*(P-1)];
	rx[i+M*N*(P-1)]=rx[i+M*(N-1)+N*M*(P-1)];

	ry[i+(N-1)*M+(P-1)*N*M]=fhy[i+(N-1)*M+(P-1)*N*M]-coefyboundFR[i-1]*vhy[i+(N-1)*M+(P-1)*N*M]+rescoefpp1boundFR[i-1]*vhy[i+M+(P-1)*M*N]+
		rescoefpm1boundFR[i-1]*vhy[i+(N-2)*M+(P-1)*M*N]+2*rescoefpp2boundFR[i-1]*vhy[i+1+(N-1)*M+(P-1)*N*M]+2*rescoefpm2boundFR[i-1]				
		*vhy[i-1+(N-1)*M+(P-1)*M*N]+rescoefpp3boundFR[i-1]*vhy[i+(N-1)*M+M*N]+rescoefpm3boundFR[i-1]*vhy[i+(N-1)*M+(P-2)*M*N]
		+viscpp1[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+1+M+(P-1)*N*M]-vhx[i-1+M+(P-1)*N*M])+viscpm1[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i-1+(N-2)*M+(P-1)*N*M]
		-vhx[i+1+(N-2)*M+(P-1)*N*M])+viscpp3[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhz[i+1+(N-1)*M+M*N]-vhz[i-1+(N-1)*M+M*N])+viscpm3[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]
		*(vhz[i-1+(N-1)*M+(P-2)*M*N]-vhz[i+1+(N-1)*M+(P-2)*M*N]);

	ry[i]=ry[i+M*(N-1)+N*M*(P-1)];
	ry[i+M*(N-1)]=ry[i+M*(N-1)+N*M*(P-1)];
	ry[i+M*N*(P-1)]=ry[i+M*(N-1)+N*M*(P-1)];
	
		
	rz[i+(N-1)*M+(P-1)*N*M]=fhz[i+(N-1)*M+(P-1)*N*M]-coefzboundFR[i-1]*vhz[i+(N-1)*M+(P-1)*N*M]+rescoefpp1boundFR[i-1]*vhz[i+M+(P-1)*M*N]+
		rescoefpm1boundFR[i-1]*vhz[i+(N-2)*M+(P-1)*M*N]+rescoefpp2boundFR[i-1]*vhz[i+1+(N-1)*M+(P-1)*N*M]+rescoefpm2boundFR[i-1]
		*vhz[i-1+(N-1)*M+(P-1)*M*N]+2*rescoefpp3boundFR[i-1]*vhz[i+(N-1)*M+M*N]+2*rescoefpm3boundFR[i-1]*vhz[i+(N-1)*M+(P-2)*M*N]
		+viscpp1[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+M+M*N]-vhx[i+M+(P-2)*N*M])+viscpm1[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhx[i+(N-2)*M+(P-2)*N*M]
		-vhx[i+(N-2)*M+M*N])+viscpp2[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]*(vhy[i+1+(N-1)*M+M*N]-vhy[i+1+(N-1)*M+(P-2)*M*N])+viscpm2[i-1+(N-2)*(M-2)+(P-2)*(M-2)*(N-1)]
		*(vhy[i-1+(N-1)*M+(P-2)*M*N]-vhy[i-1+(N-1)*M+M*N]);

	rz[i]=rz[i+M*(N-1)+N*M*(P-1)];
	rz[i+M*(N-1)]=rz[i+M*(N-1)+N*M*(P-1)];
	rz[i+M*N*(P-1)]=rz[i+M*(N-1)+N*M*(P-1)];

	}
	



plhs[0]=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS, mxREAL);
plhs[1]=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS, mxREAL);
plhs[2]=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS, mxREAL);

memcpy(mxGetPr(plhs[0]),rx,sizeof(double)*M*N*P);
memcpy(mxGetPr(plhs[1]),ry,sizeof(double)*M*N*P);
memcpy(mxGetPr(plhs[2]),rz,sizeof(double)*M*N*P);

return;
}
