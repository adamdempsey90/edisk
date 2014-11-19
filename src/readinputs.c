#include "edisk.h" 

void read_inputs(char *inputdir) {
	FILE *f;
	char garbage[100], outdir[100], ifname[100];
	double  m,
			rmin,
			rmax,
			h,
			indfl,
			alpha,
			sig0,
			indsig,
			q,
			e0,
			w0,
			rs,
			ms,
			oms,
			t0,
			tau,
			endt,
			tol;
	
	int numf;
	
	size_t len = strlen(inputdir);
	if (inputdir[len-1] != '/' && len != 0) inputdir[len] = '/';
	strcpy(ifname,inputdir);
	strcat(ifname,"params.in");
	f=fopen(ifname,"r");
	if (f==NULL) printf("\n\nERROR Can't Find Input File!\n\n");
	fgets(garbage,sizeof(garbage),f);	// Input Parameters
	fgets(garbage,sizeof(garbage),f);	// Disk Parameters
	fscanf(f,"Nr = %d \n",&NR);
	fscanf(f,"m = %lg \n",&m);
	fscanf(f,"rmin = %lg \n",&rmin);
	fscanf(f,"rmax = %lg \n",&rmax);
	fscanf(f,"h0 =  %lg \n",&h);
	fscanf(f,"flare index = %lg \n", &indfl);
	fscanf(f,"alpha =  %lg \n",&alpha);
	fscanf(f,"sigma0 =  %lg \n",&sig0);
	fscanf(f,"sigma index =  %lg \n",&indsig);
	fscanf(f,"rot index =  %lg \n",&q);
	fgets(garbage,sizeof(garbage),f);	// Initial Eccentricty
	fscanf(f,"initial e = %lg \n",&e0);
	fscanf(f,"initial a.o.p = %lg \n",&w0);
	fgets(garbage,sizeof(garbage),f);	// Star Parameters
	fscanf(f,"rsoft =  %lg \n",&rs);
	fscanf(f,"Ms =  %lg \n",&ms);
	fscanf(f,"oms =  %lg \n",&oms);
	fgets(garbage,sizeof(garbage),f);	// Time Parameters
	fscanf(f,"t0 =  %lg \n",&t0);
	fscanf(f,"tau =  %lg \n",&tau);
	fscanf(f,"endt =  %lg \n",&endt);
	fscanf(f,"numf =  %d \n",&numf);
	fscanf(f,"tol =  %lg \n",&tol);
	fscanf(f,"outputdir = %s \n",outdir);
		
	fclose(f);
	
	Params->m = m;
	Params->rmin = rmin;
	Params->rmax = rmax;
	Params->h = h;
	Params->indfl = indfl;
	Params->alpha = alpha;
	Params->sig0 = sig0;
	Params->indsig = indsig;
	Params->q = q;
	Params->e0 = e0;
	Params->w0 = w0;
	Params->rs = rs;
	Params->ms = ms;
	Params->oms = oms;
	Params->t0 = t0;
	Params->tau = tau;
	Params->endt = endt;
	Params->numf = numf;
	Params->tol = tol;
	Params->dr  = (rmax - rmin) / NR;
	Params->om0 = sqrt(Params->ms);
	strcpy(Params->outdir,outdir);
	NTOT = NR+2*NG;
	
#ifndef COMPANION
	Params->oms = 0;
#endif   
   
  	mkdir(outdir,0777);
	
	MPI_Printf("# Input Parameters #\n \
		\t# Disk Parameters	#\n \
		\tNr = %d\n \
		\tm = %lg\n \
		\trmin = %lg\n \
		\trmax = %lg\n \
		\th0 = %lg\n \
		\tflare index = %lg\n \
		\talpha = %lg\n \
		\tsigma0 = %lg\n \
		\tsigma index = %lg\n \
		\trot index =  %lg\n \
		\t# Initial Eccentricity #\n \
		\tinitial e = %lg\n \
		\tinital a.o.p = %lg\n \
		\t# Star Parameters	#\n \
		\trsoft = %lg\n \
		\tMs = %lg\n \
		\toms = %lg\n \
		\t# Time Parameters	#\n \
		\tt0 = %lg\n \
		\ttau = %lg\n \
		\tendt = %lg\n \
		\tnumf = %d\n \
		\ttol = %lg\n \
		\toutputdir = %s\n",
		NR,
		Params->m,
		Params->rmin,
		Params->rmax,
		Params->h,
		Params->indfl,
		Params->alpha,
		Params->sig0,
		Params->indsig,
		Params->q,
		Params->e0,
		Params->w0,
		Params->rs,
		Params->ms,
		Params->oms,
		Params->t0,
		Params->tau,
		Params->endt,
		Params->numf,
		Params->tol,
		Params->outdir);
		
	output_params();
	
/* Send out data to rest of processors */
	
// 	int i;
// 	if (rank==0) {	
// 		for(i=0;i<np;i++) Nxproc[i] = Nx/np;
// 		Nxproc[np-1] += Nx % np;
// 	}
	

	return;


}
