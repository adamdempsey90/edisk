#include "edisk.h"


void user_ic(Mode *fld);
void init_fld(Mode *fld) {
	int i;
	double dr = Params->dr;
	double r;
	istart = NG;
	iend = NR-NG+1;
	
	fld->sig[0] = 0;
	fld->sig[NR-NG] = 0;
	for(i=0;i<NR;i++) {
		r = (Params->rmin) + (.5 + i ) * dr;
		fld->r[i] = r;
		
		Params->hor[i] = (Params->h)*pow(r,Params->indfl);

		
		bfld->sig[i] = (Params->sig0)*pow(r,Params->indsig);
		
		bfld->omk[i] = (Params->om0)*pow(r,Params->q);
		bfld->drv[i] = -(bfld->omk[i])*(2*(Params->indfl) + Params->q);
		
		Params->c2[i] = (Params->hor[i])*
							(Params->hor[i])*r*r*(bfld->omk[i])*(bfld->omk[i]);
		
		Params->indnu = 2 *(Params->indfl) + Params->q + 2;
		Params->nu[i] = (Params->alpha)*(Params->c2[i])*(Params->hor[i])*r;
		
		
		if (Params->oms != 0) {
			bfld->omk[i] +=  (Params->c2[i])*(Params->c2[i])*(Params->indsig)
							/ (2 * (Params->oms) * r*r);
		}
		bfld->drv[i] += bfld->omk[i]*( 1 + 2*(Params->q + Params->indfl));
		bfld->dlomk[i] = Params->q + (1-(Params->om0)*pow(r,Params->q)/(bfld->omk[i]))
							*(2*Params->indfl + Params->q);
	
		bfld->v[i] = r * (bfld->omk[i]);
		bfld->u[i] = 0;
		bfld->dru[i] = 0;
		fld->u[i] = 0;
		fld->v[i] = 0;
		fld->sig[i+istart] = 0;
	}
	
	user_ic(fld);
	
	return;
}

void user_ic(Mode *fld) {
/* Set initial conditions here.
	Could be setting u,v to give initial eccentricity profile
*/
	int i;
	for(i=0;i<NR;i++) {
		
		fld->u[i] += I*(bfld->v[i])*.01;
		fld->v[i] += .5*(bfld->v[i])*.01;
		fld->sig[i+istart] += 0;
	
	}
	return;
}