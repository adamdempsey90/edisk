#include "edisk.h"

void wavekillbc(Mode *fld,double dt)
{
	int i,iflag,oflag;
	double R,tau,x,dtdtau;
	double x_in;
	if (fld->lr[istart] < 0 ) x_in = (fld->lr[istart])*.8;
	else x_in = (fld->lr[istart])*1.2;
//	const double x_in = 0;
	const double x_out = (fld->lr[iend-1])*0.8;
	const double tauin = .5/(bfld->omk[istart]);
	const double tauout = .05/(bfld->omk[iend-1]);
	double complex ubc, vbc, sbc;
	
#ifdef OPENMP
        #pragma omp parallel private(i,iflag,oflag,x,R,ubc,vbc,sbc,tau,dtdtau) shared(fld)
        #pragma omp for schedule(static)
#endif	
	for(i=istart;i<iend;i++) {
		x = fld->lr[i];
		R=0;
	
#ifdef KILLOUT
		if (x > x_out) {
/* Outer Boundary */
			R = (x-x_out)/(fld->lr[iend-1] - x_out);
			ubc = u_out_bc;
			vbc = v_out_bc;
			sbc = s_out_bc;
			tau = tauout;
		}
#endif
#ifdef KILLIN
		if (x < x_in)  {
			R = (x_in - x)/(x_in - fld->lr[istart]);
			ubc = u_in_bc;
			vbc = v_in_bc;
			sbc = s_in_bc;
			tau = tauin;
		}
#endif
		R *= R;
		

		if (R>0.0) {
			tau /= R; 
			dtdtau = dt/tau;
			
//			fld->u[i] /= (1+ dtdtau);
//			fld->v[i] /= (1+dtdtau);
			fld->sig[i] /= (1+dtdtau);
			
			
			
//			fld->u[i] = (fld->u[i] + dtdtau * ubc)/(1+dtdtau );
//			fld->v[i] = (fld->v[i]+ dtdtau * vbc)/(1+dtdtau );
//			fld->sig[i] = (fld->sig[i] + dtdtau * sbc)/(1+dtdtau);
			
			
		}
	}
	return;
}