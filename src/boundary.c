#include "edisk.h"


void set_bc(Mode *fld) {
	int i;
	
	for(i=0;i<istart;i++) {
	
#ifdef ZEROBC
		fld->u[i] = fld->u[2*istart-i];
		fld->v[i] = fld->v[2*istart-i];
		fld->sig[i] = fld->sig[2*istart-i];
		fld->u[i+iend] = fld->u[iend-i-1];
		fld->v[i+iend] = fld->v[iend-i-1];
		fld->sig[i+iend] = fld->sig[iend-i-1];
#else
		fld->u[i] = fld->u[istart];
		fld->v[i] = fld->v[istart];
		fld->sig[i] = fld->sig[istart];
		fld->u[i+iend] = fld->u[iend-1];
		fld->v[i+iend] = fld->v[iend-1];
		fld->sig[i+iend] = fld->sig[iend-1];	

#endif
	}
	return;
}

void wavekillbc(Mode *fld,double dt)
{
	int i;
	double R,tau,x,dtdtau;
	const double x_in = (fld->r[istart])*1.1;
	const double x_out = (fld->r[iend-1])*0.9;
	const double tau0 = .1/(bfld->omk[istart]);
	
	for(i=istart;i<iend;i++) {
		x = fld->r[i];
		R=0;
		if (x > x_out) R = (x-x_out)/(fld->r[iend-1] - x_out);
		if (x < x_in) R = (x_in - x)/(x_in - fld->r[istart]);

		R *= R;
		tau = tau0;

		if (R>0.0) {
			tau /= R; 
			dtdtau = dt/tau;
			fld->u[i] = (fld->u[i])/(1+dtdtau );
			fld->v[i] = (fld->v[i])/(1+dtdtau );
			fld->sig[i] = (fld->sig[i])/(1+dtdtau);				

		}
	}
	return;
}