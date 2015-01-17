#include "edisk.h"


double complex *Ceccentricity;
double *eccentricity;
double *periapse;

int maximum_ecc(void);
void calc_eccentricity(double complex *u, double complex *v, double *vybar);

void init_history(void) {
	FILE *f;
	char fname[STRLEN];
	strcpy(fname,Params->outdir);
	strcat(fname,"history.dat");
	f = fopen(fname,"w");
	
	fprintf(f,"#t\tdt\txs\tys\trs\tps\t<Ke>\t<e>\t<w>\temax\twmax\temid\n");
	

	
	fclose(f);
	
	eccentricity = (double *)malloc(sizeof(double)*NTOT);
	periapse = (double *)malloc(sizeof(double)*NTOT);
	Ceccentricity = (double complex *)malloc(sizeof(double complex)*NTOT);
	
	
	
	return;
}



void free_history(void) {

	free(eccentricity);
	free(Ceccentricity);
	free(periapse);
	
	return;
}
void history(double t, double dt, Mode *fld) {
	FILE *f;
	char fname[STRLEN];
	strcpy(fname,Params->outdir);
	strcat(fname,"history.dat");
	
	int max_ind;
	double kinetic,ebar,wbar,emax,wmax,emid,u2bar,v2bar;
	double rs,ps,xs,ys;
	double area_element = 1.0/(fld->r[iend-1]*fld->r[iend-1] - fld->r[istart]*fld->r[istart]); 
	
#ifdef INDIRECT
	rs = CentralStar->r;
	ps = CentralStar->phi;
	xs = CentralStar->x;
	ys = CentralStar->y;
#else
	rs = 0;
	ps = 0;
	xs = 0;
	ys = 0;
#endif	



	calc_eccentricity(fld->u,fld->v,bfld->v);
	
	emid = eccentricity[istart+NR/2];
	
	max_ind =  maximum_ecc();
	
	emax = eccentricity[max_ind];
	wmax = periapse[max_ind];
	
	
	
		
	u2bar = c_area_integrate(fld->r,fld->u,fld->u,istart,iend)*area_element;
	
	v2bar = c_area_integrate(fld->r,fld->v,fld->v,istart,iend)*area_element;


	ebar = r_area_integrate(fld->r,eccentricity,NULL,istart,iend)*area_element;
	wbar = r_area_integrate(fld->r,periapse,NULL,istart,iend)*area_element;
	
	
	kinetic= .5*(u2bar + v2bar);
				
	
	
	
	
	f = fopen(fname,"a");
	fprintf(f,"%lg\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
		t,dt,xs,ys,rs,ps,kinetic,ebar,wbar,emax,wmax,emid);

	fclose(f);
	
	return;
}

int maximum_ecc(void) {
	int i,j;
	
	j=0;
	for(i=istart;i<iend;i++) {
		if (eccentricity[i] > eccentricity[j]) j = i;
	}

	return j;
}

void calc_eccentricity(double complex *u, double complex *v, double *vybar) {
	int i;
	
	for(i=0;i<NTOT;i++) {
		
		Ceccentricity[i] = (2*v[i] - I*u[i])/(2*vybar[i]);
		
		eccentricity[i] = sqrt(creal(Ceccentricity[i] * conj(Ceccentricity[i])));
		
		periapse[i] = atan2(cimag(Ceccentricity[i]),creal(Ceccentricity[i]));
	}	

	return;
}

