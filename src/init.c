#include "edisk.h"

void user_ic(Mode *fld);
void init_fld(Mode *fld) {
	int i;
	double dr = Params->dr;
	double r;
	

	istart = NG;
	iend = NR+NG;
	
	fld->sig[0] = 0;
	fld->sig[NR-NG] = 0;
	

	for(i=0;i<NTOT;i++) {
		r = (Params->rmin) + (.5 + i -NG ) * dr;
		fld->r[i] = r;
		
	
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
	double e0 = .3;
	double complex w0 = 0;
	double complex E0 = e0 * cexp(I*w0);
	for(i=0;i<NTOT;i++) {
		
		fld->u[i] += I*(bfld->v[i])*E0;
		fld->v[i] += .5*(bfld->v[i])*E0;	
	}
	return;
}