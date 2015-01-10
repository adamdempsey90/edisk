#include "edisk.h"

void transport_step(double dt, Mode *fld) {
	int i;
	double complex fac;
	for(i=istart; i<iend;i++) {
//		fac = (1 + I*(fld->m)*(bfld->omk[i])*dt - .5 * (fld->m)*(bfld->omk[i])*(fld->m)*(bfld->omk[i])*dt*dt
//				- (I/6.) *pow((fld->m)*(bfld->omk[i])*dt,3));
		fac = cexp(I*(fld->m)*(bfld->omk[i])*dt);
		fld->u[i] *= fac;
		fld->v[i] *= fac;
		fld->sig[i] *= fac;
	}
	return;
}