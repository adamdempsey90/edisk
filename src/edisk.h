#include "defines.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <sys/stat.h>


#define MPI_Printf	printf

typedef struct Mode {

	double complex u,v,sig;
	double m;


} Mode;



double r0, q, om0, omf;

