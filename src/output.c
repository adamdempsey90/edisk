#include "planetdisk.h"

#define STRLEN 100
/* Option to write output in real space or complex space */

void write_header(FILE *f);
void write_cheader(FILE *f);



void output(Field *fld) {
	FILE *f;
	char fname[STRLEN];
	
	
	strcpy(fname,Params->outdir); 

	sprintf(name,"output_%d.dat",outnum);

	strcat(fname,name); 
	

// 	sprintf(fnameu,"outputs/id%d/vx_%d.dat",rank,outnum);
// 	sprintf(fnamev,"outputs/id%d/vy_%d.dat",rank,outnum);
// 	sprintf(fnames,"outputs/id%d/dens_%d.dat",rank,outnum);


#ifdef OUTBINARY
	f = fopen(fname,"wb");
	if (f == NULL) printf("ERROR: Couldn't open output file\n");

	write_cheader(f);
	
	fwrite((double *)&fld->r[0],sizeof(double),NR,f);
	fwrite((double *)&fld->u[0],sizeof(double),NR,f);
	fwrite((double *)&fld->v[0],sizeof(double),NR,f);
	fwrite((double *)&fld->sig[istart],sizeof(double),NR,f);
	fwrite((double *)&bfld->v[0],sizeof(double),NR,f);
	fwrite((double *)&bfld->sig[0],sizeof(double),NR,f);

#else
	f = fopen(fname,"w");
	if (f == NULL) printf("ERROR: Couldn't open output file\n");

	for(i=0;i<NR;i++) {
		fprintf(f,"%lg\t%lg\%lg\t%lg\t%lg\%lg\n",fld->r[i],fld->u[i],fld->v[i],
			fld->sig[i+istart],bfld->v[i],bfld->sig[i]);
	}

#endif

	fclose(f);	
	
	outnum++;
	
	return;

}


void output_params(Field *fld) {
	int i,Nxtot;
	char fname[STRLEN];
	
	for(i=0, Nxtot=0;i<np;i++) Nxtot+=Nxproc[i];
	
	strcpy(fname,fld->Params->outdir);
	strcat(fname,"params.txt");
	FILE *f = fopen(fname,"a");
	fprintf(f,"Input Parameters: \n \
		Nx = %d\n \
		Ny = %d\n \
		Lx = %lg\n \
		Ly =  %lg\n \
		dx = %lg\n \
		xsoft =  %lg\n \
		h =  %lg\n \
		Mp =  %lg\n \
		nu =  %lg\n \
		q =  %lg\n \
		omega =  %lg\n \
		sig0 = %lg\n \
		Time Parameters \n \
		t0 =  %lg\n \
		tau =  %lg\n \
		endt =  %lg\n \
		numf =  %d\n \
		tol =  %lg\n",
			  Nxtot, fld->Params->Ny, fld->Params->Lx, fld->Params->Ly, fld->Params->dx, fld->Params->xs, fld->Params->h, 
			  fld->Params->Mp, fld->Params->nu, fld->Params->q, fld->Params->omega, fld->Params->sig0,
			  fld->Params->t0, fld->Params->tau,   fld->Params->endt, fld->Params->numf,fld->Params->tol);
	fclose(f);


	return;

}
void write_header(FILE *f) {


	double dNr = (double)NR;
	double num = (double)6;
	fwrite(&dNr,sizeof(double),1,f);
	fwrite(&num,sizeof(double),1,f);

	return;
}
void init_output(char *dir) {
	char idstr[50];	
	size_t len = strlen(dir);
	if (dir[len-1] != '/') dir[len] = '/';
	mkdir(dir,0777);
	sprintf(idstr,"id%d/",rank);
	strcat(dir,idstr);
	mkdir(dir,0777);

	return;
}