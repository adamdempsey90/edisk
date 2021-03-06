#include "edisk.h" 

#ifdef LOG10
#define UNLOG(a) pow(10,a)
#else
#define UNLOG(a) exp(a)
#endif


void read_inputs(char *inputdir) {
	FILE *f;
	char garbage[100], outdir[100], ifname[100];
	char *gchar;
	double  m,
			rmin,
			rmax,
			cfl,
			h,
			indfl,
			alpha_s,
			alpha_b,
			Mdisk,
			sig0,
			mdisk,
			indsig,
			q,
			eps_sg,
			e0,
			w0,
			rs,
			init_star_rad,
			init_star_phi,
			ms,
			oms,
			t0,
			tau,
			endt,
			tol;
	
	int numf,hist_dt,read_result;
	
	size_t len = strlen(inputdir);
	if (inputdir[len-1] != '/' && len != 0) inputdir[len] = '/';
	strcpy(ifname,inputdir);
	strcat(ifname,"params.in");
	f=fopen(ifname,"r");
	if (f==NULL) printf("\n\nERROR Can't Find Input File!\n\n");
	gchar=fgets(garbage,sizeof(garbage),f);	// Input Parameters
	gchar=fgets(garbage,sizeof(garbage),f);	// Disk Parameters
	read_result=fscanf(f,"Nr = %d \n",&NR);
	read_result=fscanf(f,"m = %lg \n",&m);
	read_result=fscanf(f,"rmin = %lg \n",&rmin);
	read_result=fscanf(f,"rmax = %lg \n",&rmax);
	read_result=fscanf(f,"cfl = %lg \n",&cfl);
	read_result=fscanf(f,"h0 =  %lg \n",&h);
	read_result=fscanf(f,"flare index = %lg \n", &indfl);
	read_result=fscanf(f,"alpha shear =  %lg \n",&alpha_s);
	read_result=fscanf(f,"alpha bulk =  %lg \n",&alpha_b);
#ifdef INPUTDISKMASS
	read_result=fscanf(f,"Mdisk =  %lg \n",&mdisk);
#else
	read_result=fscanf(f,"sigma0 =  %lg \n",&sig0);
#endif
	read_result=fscanf(f,"sigma index =  %lg \n",&indsig);
	read_result=fscanf(f,"rot index =  %lg \n",&q);
	read_result=fscanf(f,"self grav soft =  %lg \n",&eps_sg);
	gchar=fgets(garbage,sizeof(garbage),f);	// Initial Eccentricty
	read_result=fscanf(f,"initial e = %lg \n",&e0);
	read_result=fscanf(f,"initial a.o.p = %lg \n",&w0);
	gchar=fgets(garbage,sizeof(garbage),f);	// Star Parameters
	read_result=fscanf(f,"rsoft =  %lg \n",&rs);
	read_result=fscanf(f,"Ms =  %lg \n",&ms);
	read_result=fscanf(f,"initial rad =  %lg \n",&init_star_rad);
	read_result=fscanf(f,"initial phi =  %lg \n",&init_star_phi);
	gchar=fgets(garbage,sizeof(garbage),f);	// Time Parameters
	read_result=fscanf(f,"t0 =  %lg \n",&t0);
	read_result=fscanf(f,"tau =  %lg \n",&tau);
	read_result=fscanf(f,"endt =  %lg \n",&endt);
	read_result=fscanf(f,"numf =  %d \n",&numf);
	read_result=fscanf(f,"hist dt =  %d \n",&hist_dt);
	read_result=fscanf(f,"tol =  %lg \n",&tol);
	read_result=fscanf(f,"outputdir = %s \n",outdir);
		
	fclose(f);
	
	Params->m = m;
	Params->rmin = rmin;
	Params->rmax = rmax;
	Params->cfl = cfl;
	Params->h = h;
	Params->indfl = indfl;
	Params->alpha_s = alpha_s;
	Params->alpha_b = alpha_b;
	Params->indsig = indsig;
	Params->q = q;
	Params->eps_sg = eps_sg;
	Params->e0 = e0;
	Params->w0 = w0;
	Params->rs = rs;
	Params->ms = ms;
	Params->init_star_rad = init_star_rad;
	Params->init_star_phi = init_star_phi;
	Params->t0 = t0;
	Params->tau = tau;
	Params->endt = endt;
	Params->numf = numf;
	Params->hist_dt = hist_dt;
	Params->tol = tol;
	Params->dr  = (rmax - rmin) / NR;
	Params->rs *= (Params->h);
	
#ifdef INPUTDISKMASS
	if (indsig != -2) {
		sig0 = mdisk * (indsig + 2) / (pow(UNLOG(rmax),indsig+2)-pow(UNLOG(rmin),indsig+2))
				/ (2*M_PI);
	}
	else {
		sig0 = mdisk / ( rmax - rmin) / (2*M_PI);
	}
#else
	if (indsig != -2) {
		mdisk = 2*M_PI*sig0 * (pow(UNLOG(rmax),indsig+2)-pow(UNLOG(rmin),indsig+2))/(indsig+2);
	}
	else {
		mdisk = 2*M_PI*sig0 * (rmax-rmin);
	}
#endif
	Params->Mdisk = mdisk;
	Params->sig0 = sig0;

	Params->om0 = sqrt(Params->ms + Params->Mdisk);
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
		\tcfl = %lg\n \
		\th0 = %lg\n \
		\tflare index = %lg\n \
		\talpha shear = %lg\n \
		\talpha bulk = %lg \n \
		\tsigma0 = %lg\n \
		\tsigma index = %lg\n \
		\tMdisk = %lg\n \
		\trot index =  %lg\n \
		\tself grav soft =  %lg\n \
		\t# Initial Eccentricity #\n \
		\tinitial e = %lg\n \
		\tinital a.o.p = %lg\n \
		\t# Star Parameters	#\n \
		\trsoft = %lg\n \
		\tMs = %lg\n \
		\tinitial rad = %lg\n \
		\tinitial phi = %lg\n \
		\t# Time Parameters	#\n \
		\tt0 = %lg\n \
		\ttau = %lg\n \
		\tendt = %lg\n \
		\tnumf = %d\n \
		\thist dt = %d\n \
		\ttol = %lg\n \
		\toutputdir = %s\n",
		NR,
		Params->m,
		UNLOG(Params->rmin),
		UNLOG(Params->rmax),
		Params->cfl,
		Params->h,
		Params->indfl,
		Params->alpha_s,
		Params->alpha_b,
		Params->sig0,
		Params->indsig,
		Params->Mdisk,
		Params->q,
		Params->eps_sg,
		Params->e0,
		Params->w0,
		Params->rs,
		Params->ms,
		Params->init_star_rad,
		Params->init_star_phi,
		Params->t0,
		Params->tau,
		Params->endt,
		Params->numf,
		Params->hist_dt,
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
