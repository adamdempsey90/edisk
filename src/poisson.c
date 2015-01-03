#include "edisk.h"

double *Kernel, *bKernel;


void poisson(Mode *fld) {
	int i,j, indxr, indxrp, indx0, indx1;
	double complex k1,k2,k3;
	double dr = Params->dr;


#ifdef OPENMP
	 #pragma omp parallel private(i,j,indxr,indxrp,indx0,indx1,k1,k2,k3) shared(fld,bfld,Kernel) 
     #pragma omp for schedule(static)
#endif
	for(i=0;i<NR;i++) {
		indxr = i+istart;

		fld->phi_sg[i] = 0;
		for(j=0 ; j < NR-1;j++) {
			indxrp = j+istart;
			indx0 = j + i*NR;
 			indx1 = j+1+i*NR;
			
			k1 = (fld->r[indxrp])*(fld->r[indxrp])*(bfld->sig[indxrp])*(fld->sig[indxrp])*Kernel[indx0];
			
			k3 = (fld->r[indxrp+1])*(fld->r[indxrp+1])*(bfld->sig[indxrp+1])*(fld->sig[indxrp+1])*Kernel[indx1];
			
			k2 = .5*(k1+k3);
			
			
			fld->phi_sg[i] += (dr/6)*( k1 + 4*k2 + k3 );
			
		}
		
		fld->gp_sg[i] = I*(fld->m)*(fld->phi_sg[i])/(fld->r[indxr]);
	}
	
	fld->gr_sg[0] = -(fld->phi_sg[1] - fld->phi_sg[0])/(dr*fld->r[istart]);
	fld->gr_sg[NR-1] = -(fld->phi_sg[NR-1]-fld->phi_sg[NR-2])/(dr*fld->r[NR-1+istart]);
	
#ifdef OPENMP
	 #pragma omp parallel private(i) shared(fld) 
     #pragma omp for schedule(static)
#endif
	for(i=1;i<NR-1;i++) {
		fld->gr_sg[i] = -(fld->phi_sg[i+1] - fld->phi_sg[i-1])/(2*dr * fld->r[i+istart]);
	}
	
// 	for(i=0;i<NR;i++) {
// 		if (i==0 || i==NR-1) {
// 			if (i==0) {
// 				fld->gr_sg[i] = -(fld->phi_sg[i]);
// 			}
// 			else {
// 				fld->gr_sg[i-1] += fld->phi_sg[i];
// 				fld->gr_sg[i-1] /= (-2*dr*(fld->r[indxr-1]));
// 				fld->gr_sg[i] = (fld->phi_sg[i] - fld->phi_sg[i-1])/ (-dr*(fld->r[indxr]));
// 			}
// 		
// 		}
// 		else {
// 		
// 			fld->gr_sg[i] = - (fld->phi_sg[i-1]);
// 			fld->gr_sg[i-1] += fld->phi_sg[i];
// 			if (i==1) {
// 				fld->gr_sg[i-1] /= (-dr*(fld->r[indxr-1]));
// 			}
// 			else {
// 				fld->gr_sg[i-1] /= (-2*dr*(fld->r[indxr-1]));
// 		
// 			}
// 		}
// 		
// 	}
	
	
	
	
	return;			


}



double greens_function(double r, double rp, double phi, double horp, double eps) {
	double chi, ans;
	
// 	chi =.25*((rorp*rorp - 2*rorp*cos(phi))/(horp*horp) + eps*eps);
// 	
// 	ans = -sqrt(2*M_PI/(horp*horp*rp*rp)) * bessk0(chi) * exp(chi);
// 	printf("%lg \t %lg \t %lg \t %lg \t %lg\n",rorp, horp,phi,chi,bessk0(chi)*exp(chi));
	
	chi = r*r + rp*rp - 2*r*rp*cos(phi) + eps*horp*horp*rp*rp;
	
	ans = 1.0/sqrt(chi);
//	printf("%lg \t %lg \t %lg \t %lg \t %lg \t %lg\n",r, rp, horp,phi,chi,ans);
	
	
	return ans;
}

double kernel_integral(double m, double r, double rp, double horp, double eps) {
	double p,dp,ans;
	double k1,k2,k3;
	p = 0;
	dp = M_PI / 200;
	
	ans = 0;
	
	while (p < M_PI) {
	
		k1 = greens_function(r,rp,p,horp,eps);
		k1 *= cos(m*p);
		
		k2 = greens_function(r,rp,p+.5*dp,horp,eps);
		k2 *= cos(m*(p+.5*dp));
		
		
		k3 = greens_function(r,rp,p+dp,horp,eps);
		k3 *= cos(m*(p+dp));
		
		
		ans += (dp/6)*(k1+4*k2+k3);
		p += dp;
	}
	
	return -2*ans;

}

void init_poisson(double m, double *r) {
	int i,j,indxr,indxrp,indx;
	Kernel = (double *)malloc(sizeof(double)*NR*NR);
	bKernel = (double *)malloc(sizeof(double)*NR*NR);
	for(i=0;i<NR;i++) {

		for(j=0;j<NR;j++) {
			indx = j + i*NR;
			indxr= i + istart;
			indxrp = j + istart;
			
			Kernel[indx] = kernel_integral(m, r[indxr],r[indxrp],
								(Params->hor[indxrp]),Params->eps_sg);
			bKernel[indx] = kernel_integral(0, r[indxr],r[indxrp],
								(Params->hor[indxrp]),Params->eps_sg);
		}
	}
	return;
}

void free_poisson(void) {
	free(Kernel);
	free(bKernel);
	return;
}


void output_selfgrav(Mode *fld) {
	int i,j;
	char fname[STRLEN];
	strcpy(fname,Params->outdir);
	strcat(fname,"selfgrav.dat");
	FILE *f = fopen(fname,"w");
	fprintf(f,"# logr \t r \t Re(phi) \t Im(phi) \t Re(gr) \t Im(gr) \t Re(gp) \t Im(gp)\n");
	for(i=0;i<NR;i++) {
		fprintf(f,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
		fld->lr[i+istart],fld->r[i+istart], bfld->phi_sg[i],bfld->gr_sg[i]);
	}
	fclose(f);
	
	strcpy(fname,Params->outdir);
	strcat(fname,"selfgrav_kernel.dat");
	f= fopen(fname,"w");
	fprintf(f,"#r_i\tr_o\tK0\tK1\n"); 
	for(i=0;i<NR;i++) {
		for(j=0;j<NR;j++) {
			fprintf(f,"%lg\t%lg\t%lg\t%lg\t",fld->r[i+istart],fld->r[j+istart],bKernel[j+i*NR],Kernel[j+i*NR]);
		}
		fprintf(f,"\n");
	}
	fclose(f);
		
	
	return;

}
