#include "edisk.h"


#define ind3(i,j) (i + 3*j)


double complex M[3][3], U[3][3], L[3][3], K[3];
double complex *G, *H, *UK, *HG;


void cranknicholson_step(double dt, Mode *fld) {
	
	int i;
	
	set_bc(fld);
	
	solver(dt,fld);

	return;
}


void solver(double dt, Mode *fld) {
	int i;
	
	int ir,ic,ind,ind4;


/* Step 1. Forward Solve */	
	for(i=istart;i<iend;i++) {
		Gpind = (i-istart - 1)*3;
		Hpind = (i-istart - 1)*3*3;
	
/* Step 1.1 Get the A,B,C, and K block matrices 	*/
		get_matrices(dt,pow(10.,fld->r[i]), fld->m,fld->nu[i],
						fld->c[i],fld->omk[i],fld->dlomk[i],
						fld->u[i-1],fld->u[i],fld->u[i+1],
						fld->v[i-1],fld->v[i],fld->v[i+1],
						fld->sig[i-1],fld->sig[i],fld->sig[i+1]);
		
/* Step 1.2 Setup UK matrix and solve for HG matrix  */

		if (i!=istart) matvec(&L,&G[gpind],&K,-1,1);
				
		for(ir=0;ir<3;i++) {
			for(ic=0;ic<4;ic++) {
				indx4 = ir + 4*ic;
				indx = ir + 3*ic;
			
				if (ic != 3) {
					UK[indx4] = -U[i][j];	
				}
				else {
					UK[indx4] = K[i];
				}
		
			}
		}	
		
		if (i != istart) matmat(&L,&H[Hpind],&M,1,1);
		solve(&M,UK)		

/* Step 1.3 Store H and G matrices */
		for(ir=0;ir<3;i++) {
			
			for(ic=0;ic<4;ic++) {
				indx4 = ir + 4*ic;
				indx = ir + 3*ic;
			
				if (ic != 3) {
					H[Hpind + indx] = UK[indx4]	
				}
				else {
					G[Gpind + ir] = UK[indx4];
				}
		
			}
		}
		

	}
/* Step 2. Backward Subsitution */

	fld->u[iend-1] = G[ (iend-1-istart)*3 + 0];
	fld->v[iend-1] = G[ (iend-1-istart)*3 + 1];
	fld->sig[iend-1] = G[ (iend-1-istart)*3 + 2];
	for(i=iend-2;i>=istart;i--) {

		
		Gind = (i-istart)*3;
		Hind = (i-istart)*3*3;
		Gnind = (i-istart+1)*3;
		
/* Step 2.1 Get solution variables, stored in G */
		matvec(&H[Hind],&G[Gnind],&G[ind]);
/* Step 3 Copy G into fld */
		fld->u[i] = G[ (i-istart)*3 + 0];
		fld->v[i] = G[ (i-istart)*3 + 1];
		fld->sig[i] = G[ (i-istart)*3 + 2];		
		

	}

	return;

}



void get_matrices(double dt, double r, double m, double nu, double c, double omk, 
							double dlomk, double complex ul, double complex uc, 
							double complex ur, double complex vl, double complex vc,
							double complex vr, double complex sl, double complex sc, 
							double complex sr) 
{
/* Fill the A,B,C,K matrices 
	A = Lower Diagonal
	B = Main Diagonal
	C = Upper Diagonal
	K = Vector of known RHS quantities
	
*/
	double complex A[3][3], B[3][3], C[3][3], F[3];
	int i,j;
	
	double dr  = (Params->dr);		// Logarithmic spacing
	double omf = Params->oms;
	double r2 = r*r;
	double dr2 = dr * dr;
	double m2 = m*m;
	
	double gam = Params->indsig + Params->indnu ;
	
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

	B[0][0] = 0;
	B[0][1] = 0;
	B[0][2] = -c;
	
	B[1][0] = 0;
	B[1][1] = 0;
	B[1][2] = 0;
	
	B[2][0]= -1.0;
	B[2][1] = 0;
	B[2][2] = 0;

/* D matrix, viscous */	

	B[0][0] += (nu/r)*( (-2./3) + (4./3)*gam);
	B[0][1] += ( - I * m *nu/(3*r));
	B[0][2] += 0;
	
	B[1][0] += (-I*m*nu/(3*r));
	B[1][1] += nu*gam/r;
	B[1][2] += omk*dlomk*nu;
	
	B[2][0] += 0;
	B[2][1] += 0;
	B[2][2] += 0;

	
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
	

/* Construct Crank-Nicholson matrices */


	for(i=0;i<3;i++) {
		for(j=0;j<3;j++) {
			ind = i + 3*j;
			
			M[ind] = .5*dt* (A[i][j] - 2*C[i][j]/dr2);
			U[ind] = .5*dt*(C[i][j]/dr2 + B[i][j] / (2*r*dr));
			L[ind] = .5*dt*(C[i][j]/dr2 - B[i][j] / (2*r*dr));
		}
	}

	K[0] =  (1+M[0][0])*uc + (1+M[0][1])*vc + (1+M[0][2])*sc
			+ U[0][0]*ur + U[0][1]*vr + U[0][2]*sr
			+ L[0][0]*ul + L[0][1]*vl + L[0][2]*sl
			+ dt*F[0];

	K[1] =  (1+M[1][0])*uc + (1+M[1][1])*vc + (1+M[1][2])*sc
			+ U[1][0]*ur + U[1][1]*vr + U[1][2]*sr
			+ L[1][0]*ul + L[1][1]*vl + L[1][2]*sl
			+ dt*F[1];
	
	K[2] =  (1+M[2][0])*uc + (1+M[2][1])*vc + (1+M[2][2])*sc
			+ U[2][0]*ur + U[2][1]*vr + U[2][2]*sr
			+ L[2][0]*ul + L[2][1]*vl + L[2][2]*sl
			+ dt*F[2];

	for(i=0;i<9;i++) {
		M[ind] = 1- M[ind];
		U[ind] *= -1;
		L[ind] *= -1;
	}

void solver_init(void) {
	
	HG = (double complex *)malloc(sizeof(double complex)*3*4);
	UK = (double complex *)malloc(sizeof(double complex)*3*4);
	
	H = (double complex *)malloc(sizeof(double complex)*3*3*NR);
	G = (double complex *)malloc(sizeof(double complex)*3*NR);
	
	return;
}

void solver_free(void) {
	
	free(HG); free(UK);
	free(H); free(G);
	return;
}