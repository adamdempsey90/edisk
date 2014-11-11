#include "edisk.h"

void fld_2_y(Mode *fld, double complex *q) {

	memcpy(&q[0],&(fld->u[0]),sizeof(double complex)*NR);
	memcpy(&q[NR],&(fld->v[0]),sizeof(double complex)*NR);
	memcpy(&q[2*NR],&(fld->sig[istart]),sizeof(double complex)*NR);
	
	return;
}
void y_2_fld(Mode *fld, double complex *q) {

	memcpy(&(fld->u[0]),&q[0],sizeof(double complex)*NR);
	memcpy(&(fld->v[0]),&q[NR],sizeof(double complex)*NR);
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