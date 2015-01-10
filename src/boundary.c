#include "edisk.h"

void user_bc(Mode *fld);

void set_bc(Mode *fld) {
	int i;
	
// 	if (creal(fld->u[iend-1]) < 0) {
// 		fld->u[iend-1] = 0;
// 	}
// 	if (creal(fld->u[istart]) > 0) {
// 		fld->u[istart] = 0;
// 	}
	
	
	for(i=0;i<istart;i++) {
	
#ifdef ZEROBC
		fld->u[i] = fld->u[2*istart-i];
		fld->v[i] = fld->v[2*istart-i];
		fld->sig[i] = fld->sig[2*istart-i];
		fld->u[i+iend] = fld->u[iend-i-1];
		fld->v[i+iend] = fld->v[iend-i-1];
		fld->sig[i+iend] = fld->sig[iend-i-1];
#else

		
//		fld->u[i] = -abs(creal(fld->u[istart]))-I*abs(cimag(fld->u[istart]));
		fld->u[i] = -abs(fld->u[istart]);
		fld->v[i] = fld->v[istart];
		fld->sig[i] = fld->sig[istart];
//		fld->sig[i] = 0;
//		fld->u[i+iend] = abs(creal(fld->u[iend-1]))+I*abs(cimag(fld->u[iend-1]));
		fld->u[i+iend] = abs(fld->u[iend-1]);
		fld->v[i+iend] = fld->v[iend-1];
//		fld->sig[i+iend] = 0;
		fld->sig[i+iend] = fld->sig[iend-1];
			
// 		fld->u[i] = u_in_bc;
// 		fld->v[i] = v_in_bc;
// 		fld->sig[i] = s_in_bc;
// 		fld->u[i+iend] = u_out_bc;
// 		fld->v[i+iend] = v_out_bc;
// 		fld->sig[i+iend] = s_out_bc;

// 		fld->u[i] = 0;
// 		fld->v[i] = 0;
// 		fld->sig[i] = 0;
// 		fld->u[i+iend] = fld->u[iend-1];
// 		fld->v[i+iend] = fld->v[iend-1];
// 		fld->sig[i+iend] = fld->sig[iend-1];
			

#endif

	}
	user_bc(fld);

	return;
}


void user_bc(Mode *fld) {
/* User can set boundary conditions here that are different than the extrapolation b.c's 
	already set (so no need to explicitly set all b.c's if they're extrapolation).
*/
	int i;
	double divv;
	
	
//	fld->u[0] = -(fld->u[1]);
//	fld->v[0] = -(fld->v[1]);
	
	
//	fld->u[iend] = u_out_bc;
//	fld->v[iend] = v_out_bc;


//	fld->sig[iend] = s_out_bc;
//	fld->u[iend] = I*(fld->m)*(fld->v[iend-1])*2*(fld->dr)+(fld->u[iend-1])*(1-2*(fld->dr));

//	fld->u[iend] = fld->u[iend-2] - 2*(fld->r[iend-1])*(Params->dr)*(fld->v[iend-1]);
//	fld->sig[iend] = 0;

//	fld->sig[iend] = - fld->sig[iend-1];

//	fld->u[istart-1] = - (fld->u[istart]);
// 	fld->v[istart-1] = - (fld->v[istart]);
 	
// 	fld->v[istart-1] = - (fld->v[istart]);
// 	fld->sig[istart-1] = - (fld->sig[istart]);



	return;
}



