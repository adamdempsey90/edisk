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
// 		fld->u[i] = u_in_bc;
// 		fld->v[i] = v_in_bc;
// 		fld->sig[i] = s_in_bc;
// 		fld->u[i+iend] = u_out_bc;
// 		fld->v[i+iend] = v_out_bc;
// 		fld->sig[i+iend] = s_out_bc;	

#endif
	}
	return;
}

void wavekillbc(Mode *fld,double dt)
{
	int i,iflag,oflag;
	double R,tau,x,dtdtau;
	double x_in;
	if (fld->r[istart] < 0 ) x_in = (fld->r[istart])*.8;
	else x_in = (fld->r[istart])*1.2;
//	const double x_in = 0;
	const double x_out = (fld->r[iend-1])*0.8;
	const double tauin = .1/(bfld->omk[istart]);
	const double tauout = .05/(bfld->omk[iend-1]);
	double complex ubc, vbc, sbc;
	
#ifdef OPENMP
        #pragma omp parallel private(i,iflag,oflag,x,R,ubc,vbc,sbc,tau,dtdtau) shared(fld)
        #pragma omp for schedule(static)
#endif	
	for(i=istart;i<iend;i++) {
		x = fld->r[i];
		R=0;
	
//#ifdef KILLOUT
		if (x > x_out) {
/* Outer Boundary */
			R = (x-x_out)/(fld->r[iend-1] - x_out);
			ubc = u_out_bc;
			vbc = v_out_bc;
			sbc = s_out_bc;
			tau = tauout;
		}
//#endif
//#ifdef KILLIN
		if (x < x_in)  {
			R = (x_in - x)/(x_in - fld->r[istart]);
			ubc = u_in_bc;
			vbc = v_in_bc;
			sbc = s_in_bc;
			tau = tauin;
		}
//#endif
		R *= R;
		

		if (R>0.0) {
			tau /= R; 
			dtdtau = dt/tau;
			
//			fld->u[i] = (fld->u[i] + dtdtau * ubc)/(1+dtdtau );
//			fld->v[i] = (fld->v[i]+ dtdtau * vbc)/(1+dtdtau );
			fld->sig[i] = (fld->sig[i] + dtdtau * sbc)/(1+dtdtau);
			
			
		}
	}
	return;
}