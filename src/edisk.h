#include "defines.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <sys/stat.h>


#define MPI_Printf	printf
#define NG 1

typedef struct Mode {

	double complex *u,*v,*sig;
	double m;
	double *r;


} Mode;

typedef struct Bmode {
	double *u, *v, *sig;
	double *drv;
	double *omk;
} Bmode;

typedef struct Parameters {
	double  m,
			rmin,
			rmax,
			h,
			indfl,
			alpha,
			sig0,
			indsig,
			q,
			rs,
			ms,
			oms,
			t0,
			tau,
			endt,
			tol;
	double dr;
	int numf;
	char outdir[100];
	
	double *hor, *nu, *c2;

} Parameters;

int NR, istart, iend;

Bmode *bfld;
Parameters *Params;