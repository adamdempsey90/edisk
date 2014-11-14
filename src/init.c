#include "edisk.h"

void user_ic(Mode *fld);
void init_fld(Mode *fld) {
	int i;
	double dr = Params->dr;
	double r,lr;
	

	istart = NG;
	iend = NR+NG;
	
	
#ifdef OPENMP
        #pragma omp parallel private(i,r,lr) shared(fld)
        #pragma omp for schedule(static)
#endif
	for(i=0;i<NTOT;i++) {
		lr = (Params->rmin) + (.5 + i -NG ) * dr;
		fld->r[i] = lr;
		r = exp(lr);
	
		Params->hor[i] = (Params->h)*pow(r,Params->indfl);

		
		bfld->sig[i] = (Params->sig0)*pow(r,Params->indsig);
		
		bfld->omk[i] = (Params->om0)*pow(r,Params->q);
		
		Params->c2[i] = (Params->hor[i])*
							(Params->hor[i])*r*r*(bfld->omk[i])*(bfld->omk[i]);
		
		Params->nu[i] = (Params->alpha)*(Params->hor[i])*(Params->hor[i])*(bfld->omk[i])*r*r;
		
		
		bfld->omk[i] *= sqrt( 1 + (Params->hor[i])*(Params->hor[i])*(Params->indsig));
		
		bfld->dlomk[i] = (Params->q) 
				+ (1-pow((Params->om0)*pow(r,Params->q)/(bfld->omk[i]),2))*(Params->indfl);
			
		bfld->v[i] = r * (bfld->omk[i]);
		bfld->u[i] = 0;
		bfld->dru[i] = 0;
		fld->u[i] = 0;
		fld->v[i] = 0;
		fld->sig[i] = 0;
	}
	Params->indnu = 2 *(Params->indfl) + Params->q + 2;
	user_ic(fld);
	
	return;
}

void user_ic(Mode *fld) {
/* Set initial conditions here.
	Could be setting u,v to give initial eccentricity profile
*/

	int i;
	double e0 = Params->e0;
	double w = Params->w0;
	double lr, r;
	double complex E0;
	double sigma = .2;
	for(i=0;i<NTOT;i++) {
		lr = (fld->r[i]);
		r = exp(lr);
		E0 = e0*cexp(I*w); //* cexp(I*drw*lr);
		fld->u[i] += I*(bfld->v[i])*E0*exp(-lr*lr/(sigma*sigma));
		fld->v[i] += .5*(bfld->v[i])*E0*exp(-lr*lr/(sigma*sigma));	
	}

/* Set B.C */	
/* Grab the inner and outer b.c's from the initialized profile. */

	u_in_bc = fld->u[istart];
//	u_in_bc = 0;
	u_out_bc = fld->u[iend-1];
	
	v_in_bc = fld->v[istart];
//	v_in_bc = 0;
	v_out_bc = fld->v[iend-1];
	
	s_in_bc = fld->sig[istart];
//	s_in_bc = 0;
	s_out_bc = fld->sig[iend-1];
	
	
	return;
}