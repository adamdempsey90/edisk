#include "edisk.h"

void user_ic(Mode *fld);
void calc_cmax(Mode *fld);
void init_fld(Mode *fld) {
	int i;
	double dr = Params->dr;
	double r,lr;
	

	istart = NG;
	iend = NR+NG;
	
	fld->m = Params->m;
	
#ifdef OPENMP
        #pragma omp parallel private(i,r,lr) shared(fld)
        #pragma omp for schedule(static)
#endif
	for(i=0;i<NTOT;i++) {
		lr = (Params->rmin) + (.5 + i -NG ) * dr;
		fld->lr[i] = lr;
#ifdef LOG10
		fld->r[i] = pow(10,lr);
#else
		fld->r[i] = exp(lr);
#endif
		r = fld->r[i];
		Params->hor[i] = (Params->h)*pow(r,Params->indfl);

		
		bfld->sig[i] = (Params->sig0)*pow(r,Params->indsig);
		
		bfld->omk[i] = (Params->om0)*pow(r,Params->q);
		bfld->dlomk[i] = (Params->q);
		
		
		Params->c2[i] = (Params->hor[i])*
							(Params->hor[i])*r*r*(bfld->omk[i])*(bfld->omk[i]);
		
		Params->nu[i] = (Params->alpha)*(Params->hor[i])*(Params->hor[i])*(bfld->omk[i])*r*r;
		
		
		bfld->omk[i] *= sqrt( 1 + (Params->hor[i])*(Params->hor[i])*(Params->indsig));
		
 		bfld->dlomk[i] += 
 				(1-pow((Params->om0)*pow(r,Params->q)/(bfld->omk[i]),2))*(Params->indfl);
			
		bfld->v[i] = r * (bfld->omk[i]);
		bfld->u[i] = 0;
		bfld->dru[i] = 0;
		fld->u[i] = 0;
		fld->v[i] = 0;
		fld->sig[i] = 0;
	}
	Params->indnu = 2 *(Params->indfl) + Params->q + 2;

	calc_cmax(fld);
	
	user_ic(fld);

#if defined(INDIRECT) || defined(COMPANION)	
	init_star(fld);
#endif
	
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
	double sigma = .05;
	double r0 = -.2;
	double aspect = (fld->lr[iend] - fld->lr[0]);
	for(i=istart;i<iend;i++) {
		lr = fld->lr[i];
		r = fld->r[i];
//		E0 = e0*cexp(I*w); //* cexp(I*drw*lr);
		
		E0 = e0 * cexp(I*w) * cos( .5*M_PI*(fld->r[iend-1] - r)/(fld->r[iend-1]-fld->r[istart]));
//		
//		E0 = e0 * cexp(I*w) * exp(-(lr-r0)*(lr-r0)/(sigma*sigma));
//		E0 = 0;

//		E0 = e0 * cexp(I*w) * (lr - fld->lr[0]) / aspect;
		fld->u[i] += I*(bfld->v[i])*E0;
		fld->v[i] += .5*(bfld->v[i])*E0;	
//		fld->sig[i] = (fld->u[i] + (fld->u[i+1] - fld->u[i-1])/(Params->dr) 
//					- I*(fld->m)*(fld->v[i]))/(I*(Params->m)*bfld->v[i]);
		fld->sig[i] = 0;
	}

/* Set B.C */	
/* Grab the inner and outer b.c's from the initialized profile. */

//	u_in_bc = fld->u[istart];
	u_in_bc = 0;
	u_out_bc = fld->u[iend-1];
	
//	v_in_bc = fld->v[istart];
	v_in_bc = 0;
	v_out_bc = fld->v[iend-1];
	
//	s_in_bc = fld->sig[istart];
	s_in_bc = 0;
	s_out_bc = fld->sig[iend-1];
//	s_out_bc = fld->sig[iend-1];
	
	for(i=0;i<istart;i++) {
		fld->u[i] = u_in_bc;
		fld->v[i] = v_in_bc;
		fld->sig[i] = s_in_bc;
		fld->u[i+iend] = u_out_bc;
		fld->v[i+iend] = v_out_bc;
		fld->sig[i+iend] = s_out_bc;
	}
	
	return;
}


void calc_cmax(Mode *fld) {
	int i;

	Params->cmax=0;
	
	for(i=istart;i<iend;i++) {
		if (sqrt(Params->c2[i]) > Params->cmax) {
			Params->cmax = sqrt(Params->c2[i]);
			Params->rcmax = fld->r[i];
		}
	}
	return;
}