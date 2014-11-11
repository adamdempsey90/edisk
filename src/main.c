#include "edisk.h"
#include <time.h>
#include <unistd.h>

int main(void) {
	int i;
	clock_t tic, toc;

	Mode *fld = (Mode *)malloc(sizeof(Mode));

	bfld = (Bmode *)malloc(sizeof(Bmode));
	Params = (Parameters *)malloc(sizeof(Parameters));
	
	
	read_inputs();
	alloc_fld(fld);
	init(fld);
	
	init_output(Params->outdir);
	
	init_rk45();

	output(fld);

  	double	h = .1;
  	double 	t=Params->t0;
  	double dt;
  	i=1;
	int term_status=0;
	int numstep=0; double avgdt=0;
	
	while (t < Params->endt)
    {
      
		dt = t;

		clear_rhs(Mode *fld);
		
		int status = rk45_step_apply(&algogas,fld,&t,&h); 
		numstep++;
		 if (status == -1) {
			MPI_Printf("ERROR With Step...\nTerminating Run...\n");
			break;
		}
		dt = t-dt;
		avgdt += dt;

		MPI_Printf ("\t step #%d, step size = %.5e, at t=%.5e \n", numstep,dt, t);
   
#ifdef WAVEKILLBC
		wavekillbc(fld,dt);
#endif
	 
	  
		if( t >= Params->t0 + i * (Params->endt) / ((double) Params->numf)) { 
			 MPI_Printf ("\t\t OUTPUT %d, step size = %.5e, at t=%.5e \n", outnum,h,t);
		
			output(fld);
			i++;
		 }
	 
  
		term_status = check_termination();
		if (term_status==-1) {
				MPI_Printf("Detected STOP file...\n");
				MPI_Printf("Outputting final state and terminating run...\n");
				output(fld);
				remove("STOP");
				break;

		}
    }
   
    toc = clock(); 
    print_time( (double)(toc - tic) / CLOCKS_PER_SEC );}
	MPI_Printf("# steps per second: %f\n", numstep /((double)(toc - tic) / CLOCKS_PER_SEC));
	MPI_Printf("Average time step: %.2e\n", avgdt/numstep);
	return 0;
}

void clear_rhs(Mode *fld) {
	int i;
	for(i=0;i<NR;i++) {
		fld->dtu[i] = 0;
		fld->dtv[i] = 0;
		fld->dtsig[i] = 0;
	}

	return;
}
int check_termination(void) {
	if( access( "STOP", F_OK ) != -1 ) return -1;
	else return 0;
}

void print_time(double t) {
	int hr, min;	
	hr = (int)floor(t/(60.*60.)); 
	t -= hr*60*60;	
	min = (int)floor(t/60);
	t -= min*60;
	
	
	if (hr==0) {
		if (min == 0) {
			printf("Total Runtime:\t%.3lgs\n",t);
			
		}
		else {
			printf("Total Runtime:\t%dm%.3lgs\n",min,t);	
		}
	}
	else {
		printf("Total Runtime:\t%dh%dm%.3lgs\n",hr,min,t);
	}
	return;
}