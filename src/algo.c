


void algogas(Mode *fld, double t) {
	int i,is;
	double complex dru,drv,drs,d2ru,d2rv,u,v,sig;
	double vr,vph,dbar,drvr,drvph,drdbar,dlomk;
	double dr,m, r,omk,omf;
	double c2,nu;
	omf = Params->oms;



	for(i=0;i<NR;i++) {
		is = i+istart;
		m = fld->m;
		dr = fld->dr;
		r = fld->r[i];
		
		drs  = (.5/dr)*(fld->sig[is+1] - fld->sig[is-1]);
		
		if (i!=0 && i !=NR-1) {
			dru  = (.5/dr)*(fld->u[i+1] - fld->u[i-1]);
			drv  = (.5/dr)*(fld->v[i+1] - fld->v[i-1]);
			d2ru = (1./(dr*dr))*(fld->u[i+1] + fld->u[i-1] - 2*fld->u[i]); 
			d2rv = (1./(dr*dr))*(fld->v[i+1] + fld->v[i-1] - 2*fld->v[i]); 
		}
		else {
			dru  = 0;
			drv = 0;
			dr2u =  0;
			dr2v = 0;
		}
		
		
		u = fld->u[i];
		v = fld->v[i];
		sig = fld->sig[is];
	
		vr = bfld->u[i];
		vph = bfld->v[i];
		dbar = bfld->sig[i];
		drvr = bfld->dru[i];
		drvph = bfld->drv[i];
		drdbar  = bfld->drsig[i];
		omk = bfld->omk[i];
		dlomk = bfld->dlomk[i];
		c2 = Params->c2[i];
		nu = Params->nu[i];
		
		
/* 	Do u eqn first */
		fld->dtu[i] += I*m*omk*u + 2*omf*v - c2*drs;

/* Viscosity */
		fld->dtu[i] += 2*nu*dr2u + 2*nu*dru/r
					-(I*nu*m/r)*(drv - v/r - I*m*u/r)
					-(2*nu/(r*r))*(u-I*m*v)
					+nu*r*2*dru*(Params->inds + Params->indnu)
					-(I*nu*m/r)*omk*dlomk*sig;
/* Stellar Potential indirect term */


/* v eom */
		fld->dtv[i] += I*m*omk*v - u*drvph - 2*omf*u + I*m*c2*sig/r;
		
/* Viscosity */
		fld->dtv[i] += nu*(dr2v + drv/r - (v+I*m*u)/(r*r) - I*m*dru/r
					-I*m*nu*(u-I*m*v)*2/(r*r) 
					+ nu*(r*drv-v-I*m*u)*(Params->inds + Params->indnu)
					+ nu*omk*dlomk*drsig;
/* Stellar Potential indirect term */

/* sig eom */
		fld->dts[i] += I*m*omk*sig - u*r*(Params->inds) - u/r - dru + I*m*v/r;
		
		
	}
	
	
	return;

} 
