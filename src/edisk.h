//#include "defines.h"
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
	double m,dr;
	double *r;
	double complex *dtu,*dtv,*dts;

} Mode;

typedef struct Bmode {
	double *u, *v, *sig;
	double *drv,*dru;
	double *omk,*dlomk;
} Bmode;

typedef struct Parameters {
	double  m,
			rmin,
			rmax,
			h,
			indfl,
			alpha,
			indnu,
			sig0,
			indsig,
			om0,
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

typedef void (*rhsfunc)(double , double complex *, double complex *, Mode *);


void algogas(double t, double complex *y, double complex *f, Mode *fld);
void init_fld(Mode *fld);
void init_output(char *dir);
void output_params(void);
void output(Mode *fld);
void read_inputs(char *inputdir);
void y_2_fld(Mode *fld, double complex *q);
void fld_2_y(Mode *fld, double complex *q);
void f_2_fld(Mode *fld, double complex *q);
void fld_2_f(Mode *fld, double complex *q);
void alloc_fld(Mode *fld);
void free_fld(Mode *fld);
void init_rk45(void);
void free_rk45(void);
int rk45_step_apply(rhsfunc func, Mode *fld,double *t, double *h);
void output_disk(double *r);
void output_rhs(Mode *fld);

int NR, istart, iend;
int outnum;

Bmode *bfld;
Parameters *Params;