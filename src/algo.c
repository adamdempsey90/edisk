#include "edisk.h"


void algogas(double t, double complex *y, double complex *f, Mode *fld) {
	int i,is;
	double complex dru,drv,drs,d2ru,d2rv,u,v,sig;
//	double vr,vph,drvr;
	double drvph,dlomk;
	double dr,m, r,omk,omf;
	double c2,nu;
	omf = Params->oms;

	
	y_2_fld(fld,y);
	dr = Params->dr;
	for(i=0;i<NR;i++) {
		is = i+istart;
		m = fld->m;
		
		r = fld->r[i];
		drs  = (.5/dr)*(fld->sig[is+1] - fld->sig[is-1]);
		
		if (i!=0 && i !=NR-1) {
			dru  = (.5/dr)*(fld->u[i+1] - fld->u[i-1]);
			drv  = (.5/dr)*(fld->v[i+1] - fld->v[i-1]);
			d2ru = (1./(dr*dr))*(fld->u[i+1] + fld->u[i-1] - 2*fld->u[i]); 
			d2rv = (1./(dr*dr))*(fld->v[i+1] + fld->v[i-1] - 2*fld->v[i]); 
		}
		else {
			dru  = 0;
			drv = 0;
			d2ru =  0;
			d2rv = 0;
		}
		
		
		u = fld->u[i];
		v = fld->v[i];
		sig = fld->sig[is];
	
//		vr = bfld->u[i];
//		vph = bfld->v[i];
//		drvr = bfld->dru[i];
		drvph = bfld->drv[i];

		omk = bfld->omk[i];
		dlomk = bfld->dlomk[i];
		c2 = Params->c2[i];
		nu = Params->nu[i];
	
	
		
/* 	Do u eqn first */
		fld->dtu[i] += I*m*omk*u + 2*omf*v - c2*drs;

/* Viscosity */
		fld->dtu[i] += 2*nu*d2ru + 2*nu*dru/r
					-(I*nu*m/r)*(drv - v/r - I*m*u/r)
					-(2*nu/(r*r))*(u-I*m*v)
					+nu*r*2*dru*(Params->indsig + Params->indnu)
					-(I*nu*m/r)*omk*dlomk*sig;
/* Stellar Potential indirect term */


/* v eom */
		fld->dtv[i] += I*m*omk*v - u*drvph - 2*omf*u + I*m*c2*sig/r;
		
/* Viscosity */
		fld->dtv[i] += nu*(d2rv + drv/r - (v+I*m*u)/(r*r) - I*m*dru/r)
					-I*m*nu*(u-I*m*v)*2/(r*r) 
					+ nu*(r*drv-v-I*m*u)*(Params->indsig + Params->indnu)
					+ nu*omk*dlomk*drs;
/* Stellar Potential indirect term */

/* sig eom */
		fld->dts[i] += I*m*omk*sig - u*r*(Params->indsig) - u/r - dru + I*m*v/r;
		
	
	}
	
	fld_2_f(fld,f);
	return;

} 
