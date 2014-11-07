#include "edisk.h"
#include <time.h>
#include <unistd.h>

int main(void) {
	int i;

	clock_t tic, toc;
	init_output(fld->Params->outdir);
	allocate_field(fld);
	init_rk45();

	output(fld);

  	double	h = .1;
  	double 	t=t0;
  	double dt;
  	i=1;
	int term_status=0;
	int numstep=0; double avgdt=0;
	
	while (t < endt)
    {
      
		dt = t;

		int status = rk45_step_apply(fld,&t,&h); 
		numstep++;
		 if (status == -1) {
			MPI_Printf("ERROR With Step...\nTerminating Run...\n");
			break;
		}
		dt = t-dt;
		avgdt += dt;

		MPI_Printf ("\t step #%d, step size = %.5e, at t=%.5e \n", numstep,dt, t);
   

		wavekillbc(fld,dt);

	 
	  
		if( t >= fld->Params->t0 + i * (fld->Params->endt) / ((double) fld->Params->numf)) { 
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