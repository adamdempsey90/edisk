#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define indx (j+N*i)
void reigenvalues(double complex *A, double complex *evals, double complex *evecs, int N);

void main(void) {
  int N = 3;
  double complex *A = (double complex *)malloc(sizeof(double complex)*N*N);
  double complex *evals = (double complex *)malloc(sizeof(double complex)*N*N);
  double complex *evecs = (double complex *)malloc(sizeof(double complex)*N*N);
  int i,j;
  
  FILE *fevals, *fevec, *fmat;
  
 
  A[0*N + 0] =  0.20521637+0.3343664*I;
  A[0*N + 1] = 0.89660690+0.36807792*I;
  A[0*N + 2] = 0.23948462+0.41853924*I;
  
  A[1*N + 0] = 0.55206727+0.91207905*I;
  A[1*N + 1] =  0.73577824+0.09491792*I;
  A[1*N + 2] =  0.52949926+0.84084549*I;
  
  A[2*N + 0] = 0.67000489+0.97893975*I;
  A[2*N + 1] =  0.71430514+0.60589646*I;
  A[2*N + 2] =  0.53604633+0.62578987*I;
  
  
  
  fmat=fopen("matrixfile.dat","w");
  printf("\nMatrix is:\n\n");
  for(i=0;i<N;i++) {
    for(j=0;j<N;j++) {
      printf("%lg + %lg I\t",creal(A[indx]), cimag(A[indx]));
      fprintf(fmat,"%lg\t%lg\t",creal(A[indx]), cimag(A[indx]));
    }
    printf("\n");
    fprintf(fmat,"\n");
  }
  
  
  reigenvalues(A,evals,evecs,N);

  fevec=fopen("evecfile.dat","w");
  printf("\nEigenvectors are :\n\n");
  for(i=0;i<N;i++) {
    for(j=0;j<N;j++) {
      printf("%lg + %lg I\t",creal(evecs[indx]), cimag(evecs[indx]));
      fprintf(fevec,"%lg\t%lg\t",creal(evecs[indx]), cimag(evecs[indx]));
    }
    printf("\n");
    fprintf(fevec,"\n");
  }
  
  
  fevals=fopen("evalsfile.dat","w");
  printf("\nEigenvalues are :\n\n");
  for(i=0;i<N;i++) {
   
    printf("%lg + I%lg\n",creal(evals[i]), cimag(evals[i]));
    fprintf(fevals,"%lg\t%lg\n",creal(evals[i]), cimag(evals[i]));
   
  }
  
  fclose(fmat); fclose(fevals); fclose(fevec);
  
  free(A); free(evecs); free(evals);
  return;
}


void reigenvalues(double complex *A, double complex *evals, double complex *evecs, int N) 
{
  int i,j;
  char JOBVL = 'N';
  char JOBVR = 'V';
  int INFO; 
  int LDA = N;
  int LDVL = N;
  int LDVR = N;
  int LWORK = 2*N;
  double RWORK[2*N];
  double complex WORK[2*N];
  
  double complex tA[N][N];
  
  for(i=0;i<N;i++) {
    for(j=0;j<N;j++) {
      tA[j][i] = A[indx];
    }
  }
  
  
  zgeev_( &JOBVL, &JOBVR, &N, &tA, &LDA, evals, NULL, &LDVL, evecs, &LDVR, &WORK, &LWORK, &RWORK, &INFO );
 
 return;
}