#include "edisk.h"

void fld_2_y(Mode *fld, double complex *q) {

	memcpy(&q[0],&(fld->u[istart]),sizeof(double complex)*NR);
	memcpy(&q[NR],&(fld->v[istart]),sizeof(double complex)*NR);
	memcpy(&q[2*NR],&(fld->sig[istart]),sizeof(double complex)*NR);
	
	return;
}
void y_2_fld(Mode *fld, double complex *q) {

	memcpy(&(fld->u[istart]),&q[0],sizeof(double complex)*NR);
	memcpy(&(fld->v[istart]),&q[NR],sizeof(double complex)*NR);
	memcpy(&(fld->sig[istart]),&q[2*NR],sizeof(double complex)*NR);
	
	return;
}

void fld_2_f(Mode *fld, double complex *q) {

	memcpy(&q[0],&(fld->dtu[0]),sizeof(double complex)*NR);
	memcpy(&q[NR],&(fld->dtv[0]),sizeof(double complex)*NR);
	memcpy(&q[2*NR],&(fld->dts[0]),sizeof(double complex)*NR);
	
	return;
}
void f_2_fld(Mode *fld, double complex *q) {

	memcpy(&(fld->dtu[0]),&q[0],sizeof(double complex)*NR);
	memcpy(&(fld->dtv[0]),&q[NR],sizeof(double complex)*NR);
	memcpy(&(fld->dts[0]),&q[2*NR],sizeof(double complex)*NR);
	
	return;
}



void matmat(double complex *A, double complex *B, double complex *C, 
					double complex alpha, double complex beta) 
{
/* Performs \alpha * A.B + \beta * C and stores the output in C.
	A,B, and C are all matrices.
	This is essenitally a wrapper for the ZGEMM BLAS routine 
*/

	char TRANSA = 'N';
	char TRANSB = 'N';
	int m = 3;
	int n = 3;
	int k = 3;
	int LDA = 3;
	int LDB = 3;
	int LDC = 3;
		 
	
	
	zgemm_(&TRANSA, &TRANSB, &m,&n,&k,&alpha,A,&LDA,B,&LDB,&beta,C,&LDC);


	return;

}


void matvec(double complex *A, double complex *B, double complex *C, 
					double complex alpha, double complex beta) 
{
/* Performs \alpha * A.B + \beta * C and stores the output in C. 
	A is a matrix, B and C are vectors.
	This is essenitally a wrapper for the ZGEMV BLAS routine 
*/

	char TRANS = 'N';
	int m = 3;
	int n = 3;
	int LDA = 3;
	int INCX = 1;
	int INCY = 1;
		 
	
	
	zgemv_(&TRANS, &m,&n,&alpha,A,&LDA,B,&INCX,&beta,C,&INCY);


	return;

}


void matsolve(double complex *A, double complex *B) {
/* Solves the system A.x = B for x 
	x is stored in B on output
	This is essentially a wrapper for the ZGESV BLAS routine
*/
	int N = 3;
	int NRHS = 4;
	int LDA = 3;
	int IPIV[N];
	int LDB = 4;
	int INFO;
	zgesv_(&N,&NRHS,A,&LDA,&IPIV,B,&LDB,&INFO);

	return;
}	