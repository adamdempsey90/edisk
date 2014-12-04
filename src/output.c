#include "edisk.h"

#define STRLEN 100
/* Option to write output in real space or complex space */

void write_header(FILE *f);
//void write_cheader(FILE *f);

int rhsnum;

void output(Mode *fld) {
	int i;
	FILE *f;
	char fname[STRLEN],name[STRLEN];
	
	
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
	fprintf(f,"#r\tRe(u)\tIm(u)\tRe(v)\tIm(v)\tRe(s)\tIm(s)\tvybar\tomk\tsigbar\n");
	for(i=0;i<NTOT;i++) {
//		MPI_Printf("writing #%d\n",i);
//		MPI_Printf("%lg\n", fld->r[i]);
// 		MPI_Printf("%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
// 			fld->r[i],
// 			creal(fld->u[i]),cimag(fld->u[i]),
// 			creal(fld->v[i]),cimag(fld->v[i]),
// 			creal(fld->sig[i]),cimag(fld->sig[i]),
// 			bfld->v[i],bfld->omk[i],bfld->sig[i]);
			
		fprintf(f,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
			fld->r[i],
			creal(fld->u[i]),cimag(fld->u[i]),
			creal(fld->v[i]),cimag(fld->v[i]),
			creal(fld->sig[i]),cimag(fld->sig[i]),
			bfld->v[i],bfld->omk[i],bfld->sig[i]);
	}

#endif

	fclose(f);	
	
	outnum++;
	
	return;

}


void output_params(void) {
	char fname[STRLEN];
		
	strcpy(fname,Params->outdir);
	strcat(fname,"params.txt");
	FILE *f = fopen(fname,"w");
	
	fprintf(f,"# Input Parameters #\n \
		\t# Disk Parameters	#\n \
		\tNr = %d\n \
		\tm = %lg\n \
		\trmin = %lg\n \
		\trmax = %lg\n \
		\tcfl = %lg\n \
		\th0 = %lg\n \
		\tflare index = %lg\n \
		\talpha = %lg\n \
		\tsigma0 = %lg\n \
		\tsigma index = %lg\n \
		\tMdisk = %lg\n \
		\trot index =  %lg\n \
		\t# Initial Eccentricity #\n \
		\tinitial e = %lg\n \
		\tinital a.o.p = %lg\n \
		\t# Star Parameters	#\n \
		\trsoft = %lg\n \
		\tMs = %lg\n \
		\t# Time Parameters	#\n \
		\tt0 = %lg\n \
		\ttau = %lg\n \
		\tendt = %lg\n \
		\tnumf = %d\n \
		\ttol = %lg\n \
		\toutputdir = %s\n",
		NR,
		Params->m,
		pow(10,Params->rmin),
		pow(10,Params->rmax),
		Params->cfl,
		Params->h,
		Params->indfl,
		Params->alpha,
		Params->sig0,
		Params->indsig,
		Params->Mdisk,
		Params->q,
		Params->e0,
		Params->w0,
		Params->rs,
		Params->ms,
		Params->t0,
		Params->tau,
		Params->endt,
		Params->numf,
		Params->tol,
		Params->outdir);
	fclose(f);


	return;

}
void output_disk(double *r) {
	int i;
	char fname[STRLEN];
	strcpy(fname,Params->outdir);
	strcat(fname,"disk.dat");
	FILE *f = fopen(fname,"w");
	fprintf(f,"# r \t h/r \t c^2 \t nu\n");
	for(i=0;i<NTOT;i++) {
		fprintf(f,"%lg\t%lg\t%lg\t%lg\n",
		r[i],
		Params->hor[i],
		Params->c2[i],
		Params->nu[i]);
	}
	fclose(f);
	return;

}

void output_rhs(Mode *fld) {
	int i;
	FILE *f;
	char fname[STRLEN],name[STRLEN];
	
	
	strcpy(fname,Params->outdir); 

	sprintf(name,"rhs_%d.dat",rhsnum);

	strcat(fname,name); 
	f = fopen(fname,"w");
	for(i=0;i<NR;i++) {
		fprintf(f,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
		fld->r[i+istart],
		creal(fld->dtu[i]),cimag(fld->dtu[i]),
		creal(fld->dtv[i]),cimag(fld->dtv[i]),
		creal(fld->dts[i]),cimag(fld->dts[i]));
	}
	rhsnum++;
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
//	char idstr[50];	
	outnum = 0; rhsnum=0;
	size_t len = strlen(dir);
	if (dir[len-1] != '/') dir[len] = '/';
	mkdir(dir,0777);
// 	sprintf(idstr,"id%d/",rank);
// 	strcat(dir,idstr);
// 	mkdir(dir,0777);

	return;
}