#include "edisk.h"

void malloc_err(const char *err_str);

void alloc_fld(Mode *fld) {
	
	fld->u = (double complex *)malloc(sizeof(double complex)*NR);
	if (fld->u == NULL) malloc_err("u");
	
	fld->v = (double complex *)malloc(sizeof(double complex)*NR);
	if (fld->v == NULL) malloc_err("v");
	
	fld->sig = (double complex *)malloc(sizeof(double complex)*(NR+2*NG));
	if (fld->sig == NULL) malloc_err("sigma");
	
	fld->dtu = (double complex *)malloc(sizeof(double complex)*NR);
	if (fld->dtu == NULL) malloc_err("dtu");
	
	fld->dtv = (double complex *)malloc(sizeof(double complex)*NR);
	if (fld->dtv == NULL) malloc_err("dtv");
	
	fld->dts = (double complex *)malloc(sizeof(double complex)*NR);
	if (fld->dts == NULL) malloc_err("dtsigma");
	
	
	fld->r = (double *)malloc(sizeof(double)*NR);
	if (fld->r == NULL) malloc_err("r");
	
	bfld->u = (double *)malloc(sizeof(double)*NR);
	if (bfld->u == NULL) malloc_err("vxbar");	

	bfld->v = (double *)malloc(sizeof(double)*NR);
	if (bfld->v == NULL) malloc_err("vybar");
	
	bfld->sig = (double *)malloc(sizeof(double)*NR);
	if (bfld->sig == NULL) malloc_err("dbar");
	
	bfld->drv = (double *)malloc(sizeof(double)*NR);
	if (bfld->drv == NULL) malloc_err("drvybar");
	
	bfld->dru = (double *)malloc(sizeof(double)*NR);
	if (bfld->dru == NULL) malloc_err("drvxbar");
	
	bfld->omk = (double *)malloc(sizeof(double)*NR);
	if (bfld->omk == NULL) malloc_err("omk");
	
	bfld->dlomk = (double *)malloc(sizeof(double)*NR);
	if (bfld->dlomk == NULL) malloc_err("dlomk");
		
	Params->hor = (double *)malloc(sizeof(double)*NR);
	if (Params->hor == NULL) malloc_err("H/R");
		
	Params->nu = (double *)malloc(sizeof(double)*NR);
	if (Params->nu == NULL) malloc_err("nu");

	Params->c2 = (double *)malloc(sizeof(double)*NR);
	if (Params->c2 == NULL) malloc_err("c2");
	return;
}


void free_fld(Mode *fld) {

	free(fld->u);
	free(fld->v);
	free(fld->sig);
	
	free(fld->r);
	
	free(fld);
	
	free(bfld->u);
	free(bfld->v);
	free(bfld->sig);

	free(bfld->drv);
	free(bfld->dru);
	free(bfld->omk);
	free(bfld->dlomk);
	
	free(bfld);
	
	free(Params->hor);
	free(Params->nu);
	free(Params->c2);
	free(Params);
	
	return;
}


void malloc_err(const char *err_str) {

	MPI_Printf("\n\n\n");
	MPI_Printf("Error Allocating:\n");
	MPI_Printf(err_str);
	MPI_Printf("\n\n\n");
	
	return;
}