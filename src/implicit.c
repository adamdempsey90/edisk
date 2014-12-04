#include "edisk.h"


#define ind3(i,j) (i + 3*j)


double complex M[3][3], U[3][3], L[3][3], K[3], UK[3][4];

typedef struct CNmats {

	double complex G[3];
	double complex H[3][3];

} CNmats;


void get_matrices(int indx, double dt, double r, double m, double nu, double c, double omk, 
							double dlomk, double complex ul, double complex uc, 
							double complex ur, double complex vl, double complex vc,
							double complex vr, double complex sl, double complex sc, 
							double complex sr);
							
void cn_solver_init(void);
void cn_solver_free(void);

void get_bc_matrix(double complex u, double complex v, double complex s);


CNmats *cn_mats;
							
void cranknicholson_step(double dt, double t, Mode *fld) {
	int i;
	
	int ir,ic;
	

/* Step 1. Forward Solve */	
	for(i=0;i<NTOT;i++) {
	
	
/* Step 1.1 Get the A,B,C, and K block matrices 	*/
		if (i!=0 && i!=NTOT-1) {
			get_matrices(i-istart,dt,pow(10.,fld->r[i]), fld->m,Params->nu[i],
							Params->c2[i],bfld->omk[i],bfld->dlomk[i],
							fld->u[i-1],fld->u[i],fld->u[i+1],
							fld->v[i-1],fld->v[i],fld->v[i+1],
							fld->sig[i-1],fld->sig[i],fld->sig[i+1]);
		}
		else {
		
			get_bc_matrix(fld->u[i],fld->v[i],fld->sig[i]);
		
		}
		
/* Step 1.2 Setup UK matrix and solve for HG matrix  */

		if (i!=0) matvec(&L[0][0],&(cn_mats[i-1].G[0]),&K[0],-1,1);



		for(ir=0;ir<3;ir++) {
			for(ic=0;ic<3;ic++) UK[ir][ic] = -U[ir][ic];
			UK[ir][3] = K[ir];
		}
			
		
		if (i != 0) matmat(&L[0][0],&(cn_mats[i-1].H[0][0]),&M[0][0],1,1);
		
		
		matsolve(&M[0][0],&UK[0][0]);		

/* Step 1.3 Store H and G matrices */
		
		for(ir=0;ir<3;ir++) {
			for(ic=0;ic<3;ic++) cn_mats[i].H[ir][ic] = UK[ir][ic];
			cn_mats[i].G[ir] = UK[ir][3];
		}

		

	}
/* Step 2. Backward Subsitution */

	cn_mats[NTOT-1].G[0] = fld->u[NTOT-1];
	cn_mats[NTOT-1].G[1] = fld->v[NTOT-1];
	cn_mats[NTOT-1].G[2] = fld->sig[NTOT-1];
	

		
	for(i=NTOT-2;i>0;i--) {
	
/* Step 2.1 Get solution variables, stored in G */	
		matvec(&(cn_mats[i].H[0][0]),&(cn_mats[i+1].G[0]),&(cn_mats[i].G[0]),1,1);

/* Step 3 Copy G into fld */		
		fld->u[i] = cn_mats[i].G[0];
		fld->v[i] = cn_mats[i].G[1];
		fld->sig[i] = cn_mats[i].G[2];
			
		

	}

	return;

}



void get_matrices(int indx, double dt, double r, double m, double nu,
						    double c, double omk, double dlomk, 
							double complex ul,double complex uc, double complex ur, 
							double complex vl, double complex vc, double complex vr, 
							double complex sl, double complex sc, double complex sr) 
{
/* Fill the A,B,C,K matrices 
	A = Main Diagonal
	B = Coefficient Matrix for D X
	C = Coefficient Matrix for D^2 X
	K = Vector of known RHS quantities
	
*/
	double complex A[3][3], B[3][3], C[3][3], F[3];
	int i,j;
	
	double dr  = (Params->dr);		// Logarithmic spacing
	double omf;
	double r2 = r*r;
	double dr2 = dr * dr;
	double m2 = m*m;
	
	double gam = Params->indsig + Params->indnu ;
	
#ifdef COMPANION
	omf = cstar->oms;
#else
	omf = 0;
#endif
	
/* Main Diagonal, inviscid */	
	A[0][0] = I*m*omk;
	A[0][1]= 2*(omf + omk);
	A[0][2] = 0;
	
	A[1][0] = -(2*omf + omk*(2+dlomk));
	A[1][1] = I*m*omk;
	A[1][2] = I*m*c/r;
	
	A[2][0] = -(Params->indsig + 1.)/r;
	A[2][1] = I*m/r;
	A[2][2] = I*m*omk;


/* Main Diagonal, viscous */
	A[0][0] += (-nu*r2)*(m2 + (4./3) + (2./3)*gam);
	A[0][1] += (nu*I*m/r2)*((7./3) + (2./3)*gam);
	A[0][2] += (-nu*I*m*omk/r) * dlomk;
	
	A[1][0] += (-nu*I*m/r2)*((7./3) - gam);
	A[1][1] += (-nu/r2)*((4./3)*m2 + gam + 1);
	A[1][2] += 0;
	
	A[2][0] += 0;
	A[2][1] += 0;
	A[2][2] += 0;


/* D matrix, inviscid */	
#ifdef IMPLICIT
	B[0][0] = 0;
	B[0][1] = 0;
	B[0][2] = -c/r;
	
	B[1][0] = 0;
	B[1][1] = 0;
	B[1][2] = 0;
	
	B[2][0]= -1.0/r;
	B[2][1] = 0;
	B[2][2] = 0;

/* D matrix, viscous */	

	B[0][0] += (nu/r2)*( (-2./3) + (4./3)*gam);
	B[0][1] += ( - I * m *nu/(3*r2));
	B[0][2] += 0;
	
	B[1][0] += (-I*m*nu/(3*r2));
	B[1][1] += nu*gam/r2;
	B[1][2] += omk*dlomk*nu/r;
	
	B[2][0] += 0;
	B[2][1] += 0;
	B[2][2] += 0;
#else
	B[0][0] = 0;
	B[0][1] = 0;
	B[0][2] = 0;
	B[1][0] = 0;
	B[1][1] = 0;
	B[1][2] = 0;
	B[2][0] = 0;
	B[2][1] = 0;
	B[2][2] = 0;
#endif	
	
/* D2 matrix, viscous */	

	C[0][0] = 4*nu/3.;
	C[0][1] = 0;
	C[0][2] = 0;
	
	C[1][0] = 0;
	C[1][1] = nu;
	C[1][2] = 0;
	
	C[2][0] = 0;
	C[2][1] = 0;
	C[2][2] = 0;

/* Force Vector */

	F[0] = 0;
	F[1] = 0;
	F[2] = 0;
	
#ifdef COMPANION 
	F[0] += cstar->gr[indx];
	F[1] += cstar->gp[indx];
#endif
	

/* Construct Crank-Nicholson matrices */

	for(i=0;i<3;i++) {
		for(j=0;j<3;j++) {
			
			M[i][j] = .5*dt* (A[i][j] - 2*C[i][j]/dr2);
			U[i][j]= .5*dt*(C[i][j]/dr2 + .5*B[i][j]/dr);
			L[i][j] = .5*dt*(C[i][j]/dr2 - .5*B[i][j]/dr);
		}
	}
	K[0] =    M[0][0]*uc + M[0][1]*vc + M[0][2]*sc
			+ U[0][0]*ur + U[0][1]*vr + U[0][2]*sr
			+ L[0][0]*ul + L[0][1]*vl + L[0][2]*sl
			+ dt*F[0] + uc;

	K[1] =    M[1][0]*uc + M[1][1]*vc + M[1][2]*sc
			+ U[1][0]*ur + U[1][1]*vr + U[1][2]*sr
			+ L[1][0]*ul + L[1][1]*vl + L[1][2]*sl
			+ dt*F[1] + vc;
	
	K[2] =    M[2][0]*uc + M[2][1]*vc + M[2][2]*sc
			+ U[2][0]*ur + U[2][1]*vr + U[2][2]*sr
			+ L[2][0]*ul + L[2][1]*vl + L[2][2]*sl
			+ dt*F[2] + sc;
	
	for(i=0;i<3;i++) {
		for(j=0;j<3;j++) {
			M[i][j] *= -1;
			U[i][j] *= -1;
			L[i][j] *= -1;
			if (i==j) M[i][j] += 1;

		}
	}

	return;
}

void get_bc_matrix(double complex u, double complex v, double complex s) {
	int i,j;
	
	for(i=0;i<3;i++) {
		for(j=0;j<3;j++) {
			U[i][j] = 0;
			L[i][j] = 0;
			if (i==j)	M[i][j] = 1;
			else M[i][j] = 0;	
		}
	}
	K[0] = u;
	K[1] = v;
	K[2] = s;

	return;
}
void cn_solver_init(void) {
	
	cn_mats = (CNmats *)malloc(sizeof(CNmats)*NTOT);
	
	
	return;
}

void cn_solver_free(void) {
	
	free(cn_mats);

	return;
}