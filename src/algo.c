#include "edisk.h"

void clear_rhs(Mode *fld);
void algogas(double t,Mode *fld);

void algo(double t, double complex *y, double complex *f, Mode *fld) {
	y_2_fld(fld,y);
	clear_rhs(fld);
	
	set_bc(fld);
	
	algogas(t,fld);
	
//	poisson(t,fld);


	fld_2_f(fld,f);

	return;
}

void clear_rhs(Mode *fld) {
	int i;
	for(i=0;i<NR;i++) {
		fld->dtu[i] = 0;
		fld->dtv[i] = 0;
		fld->dts[i] = 0;
	}

	return;
}
void algogas(double t,Mode *fld) {
	int i,it;
	double complex dru,drv,drs,d2ru,d2rv,u,v,sig;
	double complex divv;
//	double vr,vph,drvr;
	double dlomk,nusinds;
	double dr,dr2,m, r,r2,omk,omf;
	double c2,nu;
	
	omf = Params->oms;

	
	dr = Params->dr;
	dr2 = dr*dr;
	m = fld->m;
#ifdef OPENMP
        #pragma omp parallel private(i,it,dru,drv,drs,d2ru,d2rv,u,v,sig,divv,r,r2,omk,c2,nu) shared(fld,dr,dr2,omf,m) 
        #pragma omp for schedule(static)
#endif	
	for(i=istart;i<iend;i++) {
		it = i-istart;
		r = exp(fld->r[i]);
		r2 = r*r;


		drs = d1c[0]*(fld->sig[i-2]) + d1c[1]*(fld->sig[i-1]) 
				+ d1c[2]*(fld->sig[i+1]) + d1c[3]*(fld->sig[i+2]);
		drs /= dr;

		dru = d1c[0]*(fld->u[i-2]) + d1c[1]*(fld->u[i-1]) 
				+ d1c[2]*(fld->u[i+1]) + d1c[3]*(fld->u[i+2]);
		dru /= dr;

		drv = d1c[0]*(fld->v[i-2]) + d1c[1]*(fld->v[i-1]) 
				+ d1c[2]*(fld->v[i+1]) + d1c[3]*(fld->v[i+2]);
		drv /= dr;


		d2rv = d2c[0]*(fld->v[i-2]) + d2c[1]*(fld->v[i-1]) 
				+ d2c[2]*(fld->v[i]) + d2c[3]*(fld->v[i+1]) + d2c[4]*(fld->v[i+2]);
		d2rv /= dr2;

		d2ru = d2c[0]*(fld->v[i-2]) + d2c[1]*(fld->v[i-1]) 
				+ d2c[2]*(fld->v[i]) + d2c[3]*(fld->v[i+1]) + d2c[4]*(fld->v[i+2]);
		d2ru /= dr2;
// 		drs  = (.5/dr)*(fld->sig[i+1] - fld->sig[i-1]);
// 		dru  = (.5/dr)*(fld->u[i+1] - fld->u[i-1]);
// 		drv  = (.5/dr)*(fld->v[i+1] - fld->v[i-1]);
// 		d2ru = (1./dr2)*(fld->u[i+1] + fld->u[i-1] - 2*fld->u[i]); 
// 		d2rv = (1./dr2)*(fld->v[i+1] + fld->v[i-1] - 2*fld->v[i]); 	
	
// 		if (i!=0 && i !=NR-1) {
// 			drs  = (.5/dr)*(fld->sig[is+1] - fld->sig[is-1]);
// 			dru  = (.5/dr)*(fld->u[i+1] - fld->u[i-1]);
// 			drv  = (.5/dr)*(fld->v[i+1] - fld->v[i-1]);
// 			d2ru = (1./(dr*dr))*(fld->u[i+1] + fld->u[i-1] - 2*fld->u[i]); 
// 			d2rv = (1./(dr*dr))*(fld->v[i+1] + fld->v[i-1] - 2*fld->v[i]); 
// 		}
// 		else {
// 			dru  = 0;
// 			drv = 0;
// 			d2ru =  0;
// 			d2rv = 0;
// 			drs = 0;
// 		}
		
		
		u = fld->u[i];
		v = fld->v[i];
		sig = fld->sig[i];
	
		divv = (u + dru - I*m*v)/r;
		
//		vr = bfld->u[i];
//		vph = bfld->v[i];
//		drvr = bfld->dru[i];

		omk = bfld->omk[i];
		dlomk = bfld->dlomk[i];
		c2 = Params->c2[i];
		nu = Params->nu[i];
	
		nusinds = Params->indsig + Params->indnu;
	
		
/* 	Do u eqn first */


		fld->dtu[it] += I*m*omk*u + 2*(omk+omf)*v - c2*drs/r;



/* Viscosity */

 		fld->dtu[it] += nu*(d2ru - m*m*u)/r2;
	
	
// 		fld->dtu[it] += 2*nu*d2ru + 2*nu*dru/r
// 					-(I*nu*m/r)*(drv - v/r - I*m*u/r)
// 					-(2*nu/(r*r))*(u-I*m*v)
// 					+nu*r*2*dru*nusinds
// 					-(I*nu*m/r)*omk*dlomk*sig
// 					-(2.*nu/(3.*r))*(divv*nusinds + dru - u/r + r*d2ru - I*m*drv+I*m*v/r);
// 					
/* Stellar Potential indirect term */


/* v eom */


		fld->dtv[it] += I*m*omk*v - u*(2*omf +omk*(dlomk +2)) + I*m*c2*sig/r;
		
		
/* Viscosity */

 		fld->dtv[it] += nu*(d2rv- m*m*v)/r2;

// 		fld->dtv[it] += nu*(d2rv + drv/r - (v+I*m*u)/(r*r) - I*m*dru/r)
// 					-I*m*nu*(u-I*m*v)*2/(r*r) 
// 					+ nu*(r*drv-v-I*m*u)*nusinds
// 					+ nu*omk*dlomk*drs
// 					-(2.*nu/(3.*r))*(-I*m*divv);
// 					
					
/* Stellar Potential indirect term */


/* sig eom */


		fld->dts[it] += I*m*omk*sig - (u/r)*(Params->indsig) - divv;
	
	}
#ifdef OUTRHS
	output_rhs(fld);
#endif
	return;

} 
