#include "edisk.h"


void starint(double q,double m,int Nph, double *grval, double *gpval);
double grintegrand(double q, double m, double phi);
double gpintegrand(double q, double m, double phi);
void output_star(double *r);

void init_star(Mode *fld) {
	int i,it;
	double r,rsoft;
	double grval,gpval;
	double indsig = Params->indsig;
	double sig0  = Params->sig0;
	double ri =fld->r[istart];
	double ro = fld->r[iend-1];
	
	if (Params->indsig != -3) {
		cstar->r = 2*M_PI*sig0 * (pow(ro,indsig+3)-pow(ri,indsig+3)) / (indsig+3); 
	}
	else {
		cstar->r = 2*M_PI*sig0 * log(ro/ri);
	}
	
	Params->oms = (Params->om0) * pow(cstar->r,-1.5);
	
	MPI_Printf("\n\nStar at radius: %lg, with rotation rate: %lg\n",cstar->r,Params->oms);
	for(i=istart;i<iend;i++) {
		it = i - istart;
		r = fld->r[i];
		rsoft = sqrt(Params->rs*Params->rs + r*r);
		
		starint(cstar->r/rsoft,fld->m,2000,&grval,&gpval);
		cstar->gr[it] = (-2*(cstar->r)	* r *pow(rsoft,-4)) * grval;		
		cstar->gp[it] = (I*(fld->m) * 2 / (rsoft*rsoft)) * gpval; 
	}
	
	output_star(fld->r);
	return;

}


void starint(double q,double m,int Nph, double *grval, double *gpval) {
	double kr1,kr2,kr3,kr4;
	double kp1,kp2,kp3,kp4;
	double phi=0;
	double h = M_PI/Nph;
	*grval = 0;
	*gpval = 0;
	
	while (phi < M_PI) {
		kr1 = grintegrand(q,m,phi);
		kp1 = gpintegrand(q,m,phi);

		kr2 = grintegrand(q,m,phi+.5*h);
		kp2 = gpintegrand(q,m,phi+.5*h);
		kr3 = grintegrand(q,m,phi+.5*h);
		kp3 = gpintegrand(q,m,phi+.5*h);
		kr4 = grintegrand(q,m,phi+h);
		kp4 = gpintegrand(q,m,phi+h);
		*grval += (h/6.) * ( kr1 + 2*kr2 + 2*kr3 + kr4);
		*gpval += (h/6.) * ( kp1 + 2*kp2 + 2*kp3 + kp4);
		phi += h;
	}
	
	return;
}


double grintegrand(double q, double m, double phi) {
	return cos(m*phi)*(cos(phi) + q) * pow(1+q*q+2*q*cos(phi),-1.5);
}
double gpintegrand(double q, double m, double phi) {
	return cos(m*phi)*pow(1+q*q+2*q*cos(phi),-.5);
}

void output_star(double *r) {
	int i;
	char fname[100];
	strcpy(fname,Params->outdir);
	strcat(fname,"star.dat");
	FILE *f = fopen(fname,"w");
	for(i=0;i<NR;i++) {
		fprintf(f,"%lg\t%lg\t%lg\t%lg\t%lg\n",
		r[i+istart],
		creal(cstar->gr[i]),
		cimag(cstar->gr[i]),
		creal(cstar->gp[i]),
		cimag(cstar->gp[i]));
	}
	fclose(f);
	return;

}