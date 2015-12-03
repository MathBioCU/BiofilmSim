#include <math.h>
#include <string.h> /*for memcpy function*/
#include "mex.h"
#include "matrix.h"
#include <omp.h>


/*to compile with open mp, use mex relaxPRESSUREprodRB3Dper2_mex.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\LDFLAGS -fopenmp"*/

void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
int i,j,k;
double *p, *f;
double *coef, *rescoefpm1, *rescoefpm2, *rescoefpm3, *rescoefpp1, *rescoefpp2, *rescoefpp3, *coefpp1, *coefpp2, *coefpp3, *coefpm1, *coefpm2, *coefpm3;
double *rescoefpm1boundT, *rescoefpp1boundT, *rescoefpm2boundT, *rescoefpp3boundT, *rescoefpm3boundT, *rescoefpp1boundB, *rescoefpm1boundB;
double *rescoefpp2boundB, *rescoefpm3boundB, *rescoefpp3boundB, *coefpp1boundT, *coefpm1boundT, *coefpp2boundB, *coefpm2boundT, *coefpp3boundT, *coefpm3boundT;
double *coefpp1boundB, *coefpm1boundB, *coefpp3boundB, *coefboundT, *coefboundB, *coefpm3boundB;
int v1;

/* Input Arguments */
p = mxGetPr(prhs[0]);
f = mxGetPr(prhs[1]);
v1 = mxGetScalar(prhs[2]);


coef = mxGetPr(mxGetField(prhs[3],0,"coef"));
coefboundB= mxGetPr(mxGetField(prhs[3],0,"coefboundB"));
coefboundT =  mxGetPr(mxGetField(prhs[3],0,"coefboundT"));
rescoefpp1 =  mxGetPr(mxGetField(prhs[3],0,"rescoefpp1"));
rescoefpp1boundB =  mxGetPr(mxGetField(prhs[3],0,"rescoefpp1boundB"));
rescoefpp1boundT =  mxGetPr(mxGetField(prhs[3],0,"rescoefpp1boundT"));
rescoefpm1 =  mxGetPr(mxGetField(prhs[3],0,"rescoefpm1"));
rescoefpm1boundB =  mxGetPr(mxGetField(prhs[3],0,"rescoefpm1boundB"));
rescoefpm1boundT = mxGetPr(mxGetField(prhs[3],0,"rescoefpm1boundT"));
rescoefpp2 =  mxGetPr(mxGetField(prhs[3],0,"rescoefpp2"));
rescoefpp2boundB =  mxGetPr(mxGetField(prhs[3],0,"rescoefpp2boundB"));
rescoefpm2 =  mxGetPr(mxGetField(prhs[3],0,"rescoefpm2"));
rescoefpm2boundT =  mxGetPr(mxGetField(prhs[3],0,"rescoefpm2boundT"));
rescoefpp3 =  mxGetPr(mxGetField(prhs[3],0,"rescoefpp3"));
rescoefpp3boundB =  mxGetPr(mxGetField(prhs[3],0,"rescoefpp3boundB"));
rescoefpp3boundT =  mxGetPr(mxGetField(prhs[3],0,"rescoefpp3boundT"));
rescoefpm3 =  mxGetPr(mxGetField(prhs[3],0,"rescoefpm3"));
rescoefpm3boundB =  mxGetPr(mxGetField(prhs[3],0,"rescoefpm3boundB"));
rescoefpm3boundT =  mxGetPr(mxGetField(prhs[3],0,"rescoefpm3boundT"));
coefpp1 =  mxGetPr(mxGetField(prhs[3],0,"coefpp1"));
coefpm1 =  mxGetPr(mxGetField(prhs[3],0,"coefpm1"));
coefpp2 =  mxGetPr(mxGetField(prhs[3],0,"coefpp2"));
coefpm2 =  mxGetPr(mxGetField(prhs[3],0,"coefpm2"));
coefpp3 =  mxGetPr(mxGetField(prhs[3],0,"coefpp3"));
coefpm3  = mxGetPr(mxGetField(prhs[3],0,"coefpm3"));
coefpp1boundB =  mxGetPr(mxGetField(prhs[3],0,"coefpp1boundB"));
coefpm1boundB =  mxGetPr(mxGetField(prhs[3],0,"coefpm1boundB"));
coefpp2boundB =  mxGetPr(mxGetField(prhs[3],0,"coefpp2boundB"));
coefpp3boundB =  mxGetPr(mxGetField(prhs[3],0,"coefpp3boundB"));
coefpm3boundB  = mxGetPr(mxGetField(prhs[3],0,"coefpm3boundB"));
coefpp1boundT =  mxGetPr(mxGetField(prhs[3],0,"coefpp1boundT"));
coefpm1boundT =  mxGetPr(mxGetField(prhs[3],0,"coefpm1boundT"));
coefpm2boundT =  mxGetPr(mxGetField(prhs[3],0,"coefpm2boundT"));
coefpp3boundT =  mxGetPr(mxGetField(prhs[3],0,"coefpp3boundT"));
coefpm3boundT  = mxGetPr(mxGetField(prhs[3],0,"coefpm3boundT"));





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
double avg;
avg=0.0;
int p2;

omp_set_num_threads(8);

for(p2=1; p2<=v1; p2++){
#pragma omp parallel
{
#pragma omp for private(i,j) 
for(k=1; k<=P-2; k=k+2){
	for(j=1; j<=N-2; j=j+2){
		for (i=1; i<=M-2; i=i+2){
			p[i+M*j+M*N*k]=-f[i+M*j+M*N*k]/coef[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]+coefpp1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*(j+1)+M*N*k]+
				coefpm1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*(j-1)+M*N*k]+coefpp2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+1+M*j+M*N*k]+
				coefpm2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i-1+M*j+M*N*k]+coefpp3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+j*M+M*N*(k+1)]+
				coefpm3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*j+M*N*(k-1)];
			}
		}
	}


#pragma omp for private(i,j)
for(k=1; k<=P-2; k=k+2){
	for(j=2; j<=N-3; j=j+2){
		for (i=2; i<=M-3; i=i+2){		
			p[i+M*j+M*N*k]=-f[i+M*j+M*N*k]/coef[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]+coefpp1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*(j+1)+M*N*k]+
				coefpm1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*(j-1)+M*N*k]+coefpp2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+1+M*j+M*N*k]+
				coefpm2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i-1+M*j+M*N*k]+coefpp3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+j*M+M*N*(k+1)]+
				coefpm3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*j+M*N*(k-1)];
			}
		}
	}

#pragma omp for private(j)
for(k=1; k<=P-2; k=k+2){
	for(j=2; j<=N-3; j=j+2){
		p[M-1+M*j+M*N*k]=-2*f[M-1+M*j+M*N*k]/coefboundT[(j-1)+(N-1)*(k-1)]+coefpp1boundT[(j-1)+(N-1)*(k-1)]*p[M-1+M*(j+1)+M*N*k]+
			coefpm1boundT[(j-1)+(N-1)*(k-1)]*p[M-1+M*(j-1)+M*N*k]+
			coefpm2boundT[(j-1)+(N-1)*(k-1)]*p[M-2+M*j+M*N*k]+coefpp3boundT[(j-1)+(N-1)*(k-1)]*p[M-1+j*M+M*N*(k+1)]+
			coefpm3boundT[(j-1)+(N-1)*(k-1)]*p[M-1+M*j+M*N*(k-1)];
		}
	}

#pragma omp for private(j)
for(k=1; k<=P-2; k=k+2){
	for(j=2; j<=N-3; j=j+2){
		p[M*j+M*N*k]=-2*f[M*j+M*N*k]/coefboundB[(j-1)+(N-1)*(k-1)]+coefpp1boundB[(j-1)+(N-1)*(k-1)]*p[M*(j+1)+M*N*k]+
			coefpm1boundB[(j-1)+(N-1)*(k-1)]*p[M*(j-1)+M*N*k]+
			coefpp2boundB[(j-1)+(N-1)*(k-1)]*p[1+M*j+M*N*k]+coefpp3boundB[(j-1)+(N-1)*(k-1)]*p[j*M+M*N*(k+1)]+
			coefpm3boundB[(j-1)+(N-1)*(k-1)]*p[M*j+M*N*(k-1)];
		}
	}

#pragma omp for
for(k=1; k<=P-2; k=k+2){
	p[M*(N-1)+M*N*k]=-2*f[M*(N-1)+M*N*k]/coefboundB[(N-2)+(N-1)*(k-1)]+coefpp1boundB[(N-2)+(N-1)*(k-1)]*p[M+M*N*k]+
		coefpm1boundB[(N-2)+(N-1)*(k-1)]*p[M*(N-2)+M*N*k]+
		coefpp2boundB[(N-2)+(N-1)*(k-1)]*p[1+M*(N-1)+M*N*k]+coefpp3boundB[(N-2)+(N-1)*(k-1)]*p[(N-1)*M+M*N*(k+1)]+
		coefpm3boundB[(N-2)+(N-1)*(k-1)]*p[M*(N-1)+M*N*(k-1)];
	p[M*N*k]=p[M*(N-1)+M*N*k];
	}

#pragma omp for
for(k=1; k<=P-2; k=k+2){
	p[M-1+M*(N-1)+M*N*k]=-2*f[M-1+M*(N-1)+M*N*k]/coefboundT[(N-2)+(N-1)*(k-1)]+coefpp1boundT[(N-2)+(N-1)*(k-1)]*p[M-1+M+M*N*k]+
		coefpm1boundT[(N-2)+(N-1)*(k-1)]*p[M-1+M*(N-2)+M*N*k]+
		coefpm2boundT[(N-2)+(N-1)*(k-1)]*p[M-2+M*(N-1)+M*N*k]+coefpp3boundT[(N-2)+(N-1)*(k-1)]*p[M-1+(N-1)*M+M*N*(k+1)]+
		coefpm3boundT[(N-2)+(N-1)*(k-1)]*p[M-1+M*(N-1)+M*N*(k-1)];
	p[M-1+M*N*k]=p[M-1+M*(N-1)+M*N*k];
	}


#pragma omp for private(i)
for(k=1; k<=P-2; k=k+2){
	for (i=2; i<=M-3; i=i+2){		
		p[i+M*(N-1)+M*N*k]=-f[i+M*(N-1)+M*N*k]/coef[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]+coefpp1[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i+M+M*N*k]+
			coefpm1[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i+M*(N-2)+M*N*k]+coefpp2[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i+1+M*(N-1)+M*N*k]+
			coefpm2[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i-1+M*(N-1)+M*N*k]+coefpp3[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i+(N-1)*M+M*N*(k+1)]+
			coefpm3[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i+M*(N-1)+M*N*(k-1)];
		p[i+M*N*k]=p[i+M*(N-1)+M*N*k];
		}
	}

/* calculate mean and subtract from p*/

/* red odd pages */
#pragma omp for private(i,j) 
for(k=2; k<=P-3; k=k+2){
	for(j=2; j<=N-3; j=j+2){
		for (i=1; i<=M-2; i=i+2){
			p[i+M*j+M*N*k]=-f[i+M*j+M*N*k]/coef[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]+coefpp1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*(j+1)+M*N*k]+
				coefpm1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*(j-1)+M*N*k]+coefpp2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+1+M*j+M*N*k]+
				coefpm2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i-1+M*j+M*N*k]+coefpp3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+j*M+M*N*(k+1)]+
				coefpm3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*j+M*N*(k-1)];
			}
		}
	}

#pragma omp for private(i) 
for(j=2; j<=N-3; j=j+2){
	for (i=1; i<=M-2; i=i+2){
		p[i+M*j+M*N*(P-1)]=-f[i+M*j+M*N*(P-1)]/coef[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]+coefpp1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i+M*(j+1)+M*N*(P-1)]+
			coefpm1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i+M*(j-1)+M*N*(P-1)]+coefpp2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i+1+M*j+M*N*(P-1)]+
			coefpm2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i-1+M*j+M*N*(P-1)]+coefpp3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i+j*M+M*N]+
			coefpm3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i+M*j+M*N*(P-2)];
		p[i+M*j]=p[i+M*j+M*N*(P-1)];
		}
	}


#pragma omp for private(i) 
for(k=2; k<=P-3; k=k+2){
	for (i=1; i<=M-2; i=i+2){
		p[i+M*(N-1)+M*N*k]=-f[i+M*(N-1)+M*N*k]/coef[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]+coefpp1[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i+M+M*N*k]+
			coefpm1[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i+M*(N-2)+M*N*k]+coefpp2[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i+1+M*(N-1)+M*N*k]+
			coefpm2[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i-1+M*(N-1)+M*N*k]+coefpp3[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i+(N-1)*M+M*N*(k+1)]+
			coefpm3[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i+M*(N-1)+M*N*(k-1)];
		p[i+M*N*k]=p[i+M*(N-1)+M*N*k];
		}
	}

#pragma omp for 
for (i=1; i<=M-2; i=i+2){
	p[i+M*(N-1)+M*N*(P-1)]=-f[i+M*(N-1)+M*N*(P-1)]/coef[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(P-2)]+coefpp1[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(P-2)]*p[i+M+M*N*(P-1)]+
		coefpm1[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(P-2)]*p[i+M*(N-2)+M*N*(P-1)]+coefpp2[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(P-2)]*p[i+1+M*(N-1)+M*N*(P-1)]+
		coefpm2[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(P-2)]*p[i-1+M*(N-1)+M*N*(P-1)]+coefpp3[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(P-2)]*p[i+(N-1)*M+M*N]+
		coefpm3[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(P-2)]*p[i+M*(N-1)+M*N*(P-2)];
	p[i+M*N*(P-1)]=p[i+M*(N-1)+M*N*(P-1)];
	p[i]=p[i+M*(N-1)+M*N*(P-1)];
	p[i+M*(N-1)]=p[i+M*N*(P-1)];
	}


#pragma omp for private(i,j) 
for(k=2; k<=P-3; k=k+2){
	for(j=1; j<=N-2; j=j+2){
		for (i=2; i<=M-3; i=i+2){
			p[i+M*j+M*N*k]=-f[i+M*j+M*N*k]/coef[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]+coefpp1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*(j+1)+M*N*k]+
				coefpm1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*(j-1)+M*N*k]+coefpp2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+1+M*j+M*N*k]+
				coefpm2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i-1+M*j+M*N*k]+coefpp3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+j*M+M*N*(k+1)]+
				coefpm3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*j+M*N*(k-1)];
			}
		}
	}

#pragma omp for private(i) 
for(j=1; j<=N-2; j=j+2){
	for (i=2; i<=M-3; i=i+2){
		p[i+M*j+M*N*(P-1)]=-f[i+M*j+M*N*(P-1)]/coef[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]+coefpp1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i+M*(j+1)+M*N*(P-1)]+
			coefpm1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i+M*(j-1)+M*N*(P-1)]+coefpp2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i+1+M*j+M*N*(P-1)]+
			coefpm2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i-1+M*j+M*N*(P-1)]+coefpp3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i+j*M+M*N]+
			coefpm3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i+M*j+M*N*(P-2)];
		p[i+M*j]=p[i+M*j+M*N*(P-1)];
		}
	}

#pragma omp for private(j)
for(k=2; k<=P-3; k=k+2){
	for(j=1; j<=N-2; j=j+2){
		p[M-1+M*j+M*N*k]=-2*f[M-1+M*j+M*N*k]/coefboundT[(j-1)+(N-1)*(k-1)]+coefpp1boundT[(j-1)+(N-1)*(k-1)]*p[M-1+M*(j+1)+M*N*k]+
			coefpm1boundT[(j-1)+(N-1)*(k-1)]*p[M-1+M*(j-1)+M*N*k]+
			coefpm2boundT[(j-1)+(N-1)*(k-1)]*p[M-2+M*j+M*N*k]+coefpp3boundT[(j-1)+(N-1)*(k-1)]*p[M-1+j*M+M*N*(k+1)]+
			coefpm3boundT[(j-1)+(N-1)*(k-1)]*p[M-1+M*j+M*N*(k-1)];
		}
	}

#pragma omp for private(j)
for(k=2; k<=P-3; k=k+2){
	for(j=1; j<=N-2; j=j+2){
		p[M*j+M*N*k]=-2*f[M*j+M*N*k]/coefboundB[(j-1)+(N-1)*(k-1)]+coefpp1boundB[(j-1)+(N-1)*(k-1)]*p[M*(j+1)+M*N*k]+
			coefpm1boundB[(j-1)+(N-1)*(k-1)]*p[M*(j-1)+M*N*k]+
			coefpp2boundB[(j-1)+(N-1)*(k-1)]*p[1+M*j+M*N*k]+coefpp3boundB[(j-1)+(N-1)*(k-1)]*p[j*M+M*N*(k+1)]+
			coefpm3boundB[(j-1)+(N-1)*(k-1)]*p[M*j+M*N*(k-1)];
		}
	}

#pragma omp for 
for(j=1; j<=N-2; j=j+2){
	p[M*j+M*N*(P-1)]=-2*f[M*j+M*N*(P-1)]/coefboundB[j-1+(N-1)*(P-2)]+coefpp1boundB[(j-1)+(N-1)*(P-2)]*p[M*(j+1)+M*N*(P-1)]+
		coefpm1boundB[(j-1)+(N-1)*(P-2)]*p[M*(j-1)+M*N*(P-1)]+
		coefpp2boundB[(j-1)+(N-1)*(P-2)]*p[1+M*j+M*N*(P-1)]+coefpp3boundB[(j-1)+(N-1)*(P-2)]*p[j*M+M*N]+
		coefpm3boundB[(j-1)+(N-1)*(P-2)]*p[M*j+M*N*(P-2)];
	p[M*j]=p[M*j+M*N*(P-1)];
	}

#pragma omp for 
for(j=1; j<=N-2; j=j+2){
	p[M-1+M*j+M*N*(P-1)]=-2*f[M-1+M*j+M*N*(P-1)]/coefboundT[j-1+(N-1)*(P-2)]+coefpp1boundT[(j-1)+(N-1)*(P-2)]*p[M-1+M*(j+1)+M*N*(P-1)]+
		coefpm1boundT[(j-1)+(N-1)*(P-2)]*p[M-1+M*(j-1)+M*N*(P-1)]+
		coefpm2boundT[(j-1)+(N-1)*(P-2)]*p[M-2+M*j+M*N*(P-1)]+coefpp3boundT[(j-1)+(N-1)*(P-2)]*p[M-1+j*M+M*N]+
		coefpm3boundT[(j-1)+(N-1)*(P-2)]*p[M-1+M*j+M*N*(P-2)];
	p[M-1+M*j]=p[M-1+M*j+M*N*(P-1)];
	}

/*black even pages*/
#pragma omp for private(i,j) 
for(k=1; k<=P-2; k=k+2){
	for(j=2; j<=N-3; j=j+2){
		for (i=1; i<=M-2; i=i+2){
			p[i+M*j+M*N*k]=-f[i+M*j+M*N*k]/coef[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]+coefpp1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*(j+1)+M*N*k]+
				coefpm1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*(j-1)+M*N*k]+coefpp2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+1+M*j+M*N*k]+
				coefpm2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i-1+M*j+M*N*k]+coefpp3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+j*M+M*N*(k+1)]+
				coefpm3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*j+M*N*(k-1)];
			}
		}
	}

#pragma omp for private(i)
for(k=1; k<=P-2; k=k+2){
	for (i=1; i<=M-2; i=i+2){		
		p[i+M*(N-1)+M*N*k]=-f[i+M*(N-1)+M*N*k]/coef[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]+coefpp1[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i+M+M*N*k]+
			coefpm1[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i+M*(N-2)+M*N*k]+coefpp2[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i+1+M*(N-1)+M*N*k]+
			coefpm2[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i-1+M*(N-1)+M*N*k]+coefpp3[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i+(N-1)*M+M*N*(k+1)]+
			coefpm3[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i+M*(N-1)+M*N*(k-1)];
		p[i+M*N*k]=p[i+M*(N-1)+M*N*k];
		}
	}

#pragma omp for private(i,j) 
for(k=1; k<=P-2; k=k+2){
	for(j=1; j<=N-2; j=j+2){
		for (i=2; i<=M-3; i=i+2){
			p[i+M*j+M*N*k]=-f[i+M*j+M*N*k]/coef[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]+coefpp1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*(j+1)+M*N*k]+
				coefpm1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*(j-1)+M*N*k]+coefpp2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+1+M*j+M*N*k]+
				coefpm2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i-1+M*j+M*N*k]+coefpp3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+j*M+M*N*(k+1)]+
				coefpm3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*j+M*N*(k-1)];
			}
		}
	}

#pragma omp for private(j)
for(k=1; k<=P-2; k=k+2){
	for(j=1; j<=N-2; j=j+2){
		p[M-1+M*j+M*N*k]=-2*f[M-1+M*j+M*N*k]/coefboundT[(j-1)+(N-1)*(k-1)]+coefpp1boundT[(j-1)+(N-1)*(k-1)]*p[M-1+M*(j+1)+M*N*k]+
			coefpm1boundT[(j-1)+(N-1)*(k-1)]*p[M-1+M*(j-1)+M*N*k]+
			coefpm2boundT[(j-1)+(N-1)*(k-1)]*p[M-2+M*j+M*N*k]+coefpp3boundT[(j-1)+(N-1)*(k-1)]*p[M-1+j*M+M*N*(k+1)]+
			coefpm3boundT[(j-1)+(N-1)*(k-1)]*p[M-1+M*j+M*N*(k-1)];
		}
	}

#pragma omp for private(j)
for(k=1; k<=P-2; k=k+2){
	for(j=1; j<=N-2; j=j+2){
		p[M*j+M*N*k]=-2*f[M*j+M*N*k]/coefboundB[(j-1)+(N-1)*(k-1)]+coefpp1boundB[(j-1)+(N-1)*(k-1)]*p[M*(j+1)+M*N*k]+
			coefpm1boundB[(j-1)+(N-1)*(k-1)]*p[M*(j-1)+M*N*k]+
			coefpp2boundB[(j-1)+(N-1)*(k-1)]*p[1+M*j+M*N*k]+coefpp3boundB[(j-1)+(N-1)*(k-1)]*p[j*M+M*N*(k+1)]+
			coefpm3boundB[(j-1)+(N-1)*(k-1)]*p[M*j+M*N*(k-1)];
		}
	}


/*black odd pages*/
#pragma omp for private(i,j) 
for(k=2; k<=P-3; k=k+2){
	for(j=2; j<=N-3; j=j+2){
		for (i=2; i<=M-3; i=i+2){
			p[i+M*j+M*N*k]=-f[i+M*j+M*N*k]/coef[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]+coefpp1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*(j+1)+M*N*k]+
				coefpm1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*(j-1)+M*N*k]+coefpp2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+1+M*j+M*N*k]+
				coefpm2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i-1+M*j+M*N*k]+coefpp3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+j*M+M*N*(k+1)]+
				coefpm3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*j+M*N*(k-1)];
			}
		}
	}

#pragma omp for private(j)
for(k=2; k<=P-3; k=k+2){
	for(j=2; j<=N-3; j=j+2){
		p[M-1+M*j+M*N*k]=-2*f[M-1+M*j+M*N*k]/coefboundT[(j-1)+(N-1)*(k-1)]+coefpp1boundT[(j-1)+(N-1)*(k-1)]*p[M-1+M*(j+1)+M*N*k]+
			coefpm1boundT[(j-1)+(N-1)*(k-1)]*p[M-1+M*(j-1)+M*N*k]+
			coefpm2boundT[(j-1)+(N-1)*(k-1)]*p[M-2+M*j+M*N*k]+coefpp3boundT[(j-1)+(N-1)*(k-1)]*p[M-1+j*M+M*N*(k+1)]+
			coefpm3boundT[(j-1)+(N-1)*(k-1)]*p[M-1+M*j+M*N*(k-1)];
		}
	}

#pragma omp for private(j)
for(k=2; k<=P-3; k=k+2){
	for(j=2; j<=N-3; j=j+2){
		p[M*j+M*N*k]=-2*f[M*j+M*N*k]/coefboundB[(j-1)+(N-1)*(k-1)]+coefpp1boundB[(j-1)+(N-1)*(k-1)]*p[M*(j+1)+M*N*k]+
			coefpm1boundB[(j-1)+(N-1)*(k-1)]*p[M*(j-1)+M*N*k]+
			coefpp2boundB[(j-1)+(N-1)*(k-1)]*p[1+M*j+M*N*k]+coefpp3boundB[(j-1)+(N-1)*(k-1)]*p[j*M+M*N*(k+1)]+
			coefpm3boundB[(j-1)+(N-1)*(k-1)]*p[M*j+M*N*(k-1)];
		}
	}

#pragma omp for
for(j=2; j<=N-3; j=j+2){
	p[M*j+M*N*(P-1)]=-2*f[M*j+M*N*(P-1)]/coefboundB[j-1+(N-1)*(P-2)]+coefpp1boundB[(j-1)+(N-1)*(P-2)]*p[M*(j+1)+M*N*(P-1)]+
		coefpm1boundB[(j-1)+(N-1)*(P-2)]*p[M*(j-1)+M*N*(P-1)]+
		coefpp2boundB[(j-1)+(N-1)*(P-2)]*p[1+M*j+M*N*(P-1)]+coefpp3boundB[(j-1)+(N-1)*(P-2)]*p[j*M+M*N]+
		coefpm3boundB[(j-1)+(N-1)*(P-2)]*p[M*j+M*N*(P-2)];
	p[M*j]=p[M*j+M*N*(P-1)];
	}

#pragma omp for
for(j=2; j<=N-3; j=j+2){
	p[M-1+M*j+M*N*(P-1)]=-2*f[M-1+M*j+M*N*(P-1)]/coefboundT[j-1+(N-1)*(P-2)]+coefpp1boundT[(j-1)+(N-1)*(P-2)]*p[M-1+M*(j+1)+M*N*(P-1)]+
		coefpm1boundT[(j-1)+(N-1)*(P-2)]*p[M-1+M*(j-1)+M*N*(P-1)]+
		coefpm2boundT[(j-1)+(N-1)*(P-2)]*p[M-2+M*j+M*N*(P-1)]+coefpp3boundT[(j-1)+(N-1)*(P-2)]*p[M-1+j*M+M*N]+
		coefpm3boundT[(j-1)+(N-1)*(P-2)]*p[M-1+M*j+M*N*(P-2)];
	p[M-1+M*j]=p[M-1+M*j+M*N*(P-1)];
	}

#pragma omp for private(i)
for(k=2; k<=P-3; k=k+2){
	for (i=2; i<=M-3; i=i+2){		
		p[i+M*(N-1)+M*N*k]=-f[i+M*(N-1)+M*N*k]/coef[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]+coefpp1[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i+M+M*N*k]+
			coefpm1[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i+M*(N-2)+M*N*k]+coefpp2[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i+1+M*(N-1)+M*N*k]+
			coefpm2[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i-1+M*(N-1)+M*N*k]+coefpp3[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i+(N-1)*M+M*N*(k+1)]+
			coefpm3[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(k-1)]*p[i+M*(N-1)+M*N*(k-1)];
		p[i+M*N*k]=p[i+M*(N-1)+M*N*k];
		}
	}
#pragma omp for
for (i=2; i<=M-3; i=i+2){
	p[i+M*(N-1)+M*N*(P-1)]=-f[i+M*(N-1)+M*N*(P-1)]/coef[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(P-2)]+coefpp1[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(P-2)]*p[i+M+M*N*(P-1)]+
		coefpm1[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(P-2)]*p[i+M*(N-2)+M*N*(P-1)]+coefpp2[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(P-2)]*p[i+1+M*(N-1)+M*N*(P-1)]+
		coefpm2[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(P-2)]*p[i-1+M*(N-1)+M*N*(P-1)]+coefpp3[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(P-2)]*p[i+(N-1)*M+M*N]+
		coefpm3[i-1+(M-2)*(N-2)+(M-2)*(N-1)*(P-2)]*p[i+M*(N-1)+M*N*(P-2)];
	p[i+M*N*(P-1)]=p[i+M*(N-1)+M*N*(P-1)];
	p[i]=p[i+M*N*(P-1)];
	p[i+M*(N-1)]=p[i+M*N*(P-1)];
	}


/*creases */
#pragma omp for
for(k=2; k<=P-3; k=k+2){
	p[M*(N-1)+M*N*k]=-2*f[M*(N-1)+M*N*k]/coefboundB[(N-2)+(N-1)*(k-1)]+coefpp1boundB[(N-2)+(N-1)*(k-1)]*p[M+M*N*k]+
		coefpm1boundB[(N-2)+(N-1)*(k-1)]*p[M*(N-2)+M*N*k]+
		coefpp2boundB[(N-2)+(N-1)*(k-1)]*p[1+M*(N-1)+M*N*k]+coefpp3boundB[(N-2)+(N-1)*(k-1)]*p[(N-1)*M+M*N*(k+1)]+
		coefpm3boundB[(N-2)+(N-1)*(k-1)]*p[M*(N-1)+M*N*(k-1)];
	p[M*N*k]=p[M*(N-1)+M*N*k];
	}
#pragma omp single
{
p[M*(N-1)+M*N*(P-1)]=-2*f[M*(N-1)+M*N*(P-1)]/coefboundB[(N-2)+(N-1)*(P-2)]+coefpp1boundB[(N-2)+(N-1)*(P-2)]*p[M+M*N*(P-1)]+
	coefpm1boundB[(N-2)+(N-1)*(P-2)]*p[M*(N-2)+M*N*(P-1)]+
	coefpp2boundB[(N-2)+(N-1)*(P-2)]*p[1+M*(N-1)+M*N*(P-1)]+coefpp3boundB[(N-2)+(N-1)*(P-2)]*p[(N-1)*M+M*N]+
	coefpm3boundB[(N-2)+(N-1)*(P-2)]*p[M*(N-1)+M*N*(P-2)];
p[M*N*(P-1)]=p[M*(N-1)+M*N*(P-1)];
p[0]=p[M*N*(P-1)];
p[M*(N-1)]=p[0];
}


#pragma omp for
for(k=2; k<=P-3; k=k+2){
	p[M-1+M*(N-1)+M*N*k]=-2*f[M-1+M*(N-1)+M*N*k]/coefboundT[(N-2)+(N-1)*(k-1)]+coefpp1boundT[(N-2)+(N-1)*(k-1)]*p[M-1+M+M*N*k]+
		coefpm1boundT[(N-2)+(N-1)*(k-1)]*p[M-1+M*(N-2)+M*N*k]+
		coefpm2boundT[(N-2)+(N-1)*(k-1)]*p[M-2+M*(N-1)+M*N*k]+coefpp3boundT[(N-2)+(N-1)*(k-1)]*p[M-1+(N-1)*M+M*N*(k+1)]+
		coefpm3boundT[(N-2)+(N-1)*(k-1)]*p[M-1+M*(N-1)+M*N*(k-1)];
	p[M-1+M*N*k]=p[M-1+M*(N-1)+M*N*k];
	}
#pragma omp single
{
p[M-1+M*(N-1)+M*N*(P-1)]=-2*f[M-1+M*(N-1)+M*N*(P-1)]/coefboundT[(N-2)+(N-1)*(P-2)]+coefpp1boundT[(N-2)+(N-1)*(P-2)]*p[M-1+M+M*N*(P-1)]+
	coefpm1boundT[(N-2)+(N-1)*(P-2)]*p[M-1+M*(N-2)+M*N*(P-1)]+
	coefpm2boundT[(N-2)+(N-1)*(P-2)]*p[M-2+M*(N-1)+M*N*(P-1)]+coefpp3boundT[(N-2)+(N-1)*(P-2)]*p[M-1+(N-1)*M+M*N]+
	coefpm3boundT[(N-2)+(N-1)*(P-2)]*p[M-1+M*(N-1)+M*N*(P-2)];
p[M-1+M*N*(P-1)]=p[M-1+M*(N-1)+M*N*(P-1)];
p[M-1]=p[M-1+M*N*(P-1)];
p[M-1+M*(N-1)]=p[M-1];
}


#pragma omp for private(i) 
for(j=2; j<=N-3; j=j+2){
	for (i=2; i<=M-3; i=i+2){
		p[i+M*j+M*N*(P-1)]=-f[i+M*j+M*N*(P-1)]/coef[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]+coefpp1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i+M*(j+1)+M*N*(P-1)]+
			coefpm1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i+M*(j-1)+M*N*(P-1)]+coefpp2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i+1+M*j+M*N*(P-1)]+
			coefpm2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i-1+M*j+M*N*(P-1)]+coefpp3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i+j*M+M*N]+
			coefpm3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i+M*j+M*N*(P-2)];
		p[i+M*j]=p[i+M*j+M*N*(P-1)];
		}
	}

#pragma omp for private(i,j) 
for(k=2; k<=P-3; k=k+2){
	for(j=1; j<=N-2; j=j+2){
		for (i=1; i<=M-2; i=i+2){
			p[i+M*j+M*N*k]=-f[i+M*j+M*N*k]/coef[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]+coefpp1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*(j+1)+M*N*k]+
				coefpm1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*(j-1)+M*N*k]+coefpp2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+1+M*j+M*N*k]+
				coefpm2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i-1+M*j+M*N*k]+coefpp3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+j*M+M*N*(k+1)]+
				coefpm3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(k-1)]*p[i+M*j+M*N*(k-1)];
			}
		}
	}


#pragma omp for private(i) 
for(j=1; j<=N-2; j=j+2){
	for (i=1; i<=M-2; i=i+2){
		p[i+M*j+M*N*(P-1)]=-f[i+M*j+M*N*(P-1)]/coef[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]+coefpp1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i+M*(j+1)+M*N*(P-1)]+
			coefpm1[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i+M*(j-1)+M*N*(P-1)]+coefpp2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i+1+M*j+M*N*(P-1)]+
			coefpm2[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i-1+M*j+M*N*(P-1)]+coefpp3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i+j*M+M*N]+
			coefpm3[i-1+(M-2)*(j-1)+(M-2)*(N-1)*(P-2)]*p[i+M*j+M*N*(P-2)];
		p[i+M*j]=p[i+M*j+M*N*(P-1)];
		}
	}

/*compute mean and subtract to keep mean(p)=i*/
/*for(k=0; k<=P-2; k++){
	for(j=0;j<=N-2; j++){
		for(i=0; i<=M-1; i++){
			avg = avg+p[i+M*j+N*M*k];
			}
		}
	}
printf("%f\n", avg);
#pragma OMP critical
{
avg=avg/((double)M*(double)N*(double)P);
}
#pragma OMP barrier

printf("%f\n",avg);
#pragma omp for private(i,j)
for(k=0; k<=P-1; k++){
	for(j=0;j<=N-1; j++){
		for(i=0; i<=M-1; i++){
		p[i+M*j+N*M*k]=p[i+M*j+N*M*k]-avg;
			}
		}
	}

*/
}
}
plhs[0]=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS, mxREAL);
memcpy(mxGetPr(plhs[0]),p,sizeof(double)*M*N*P);

return;
};
