


void algo(Modes *fld, double t) {
	int i;
	double complex dru,drv,drsig,d2ru,d2rv,d2rsig,u,v,sig;
	double dx,m, rc,omk;
	for(i=istart;i<iend;i++) {
		m = fld->m;
		dx = fld->dr;
		rc = fld->r[i];
		dru  = (.5/dr)*(fld->u[i+1] - fld->u[i-1]);
		drv  = (.5/dr)*(fld->v[i+1] - fld->v[i-1]);
		drsig  = (.5/dr)*(fld->sig[i+1] - fld->sig[i-1]);
		drpres = (.5/dr)*(fld->pres[i+1] - fld->pres[i-1]);
		d2ru = (1./dr)*(fld->u[i+1] + fld->u[i-1] - 2*fld->u[i]); 
		d2rv = (1./dr)*(fld->v[i+1] + fld->v[i-1] - 2*fld->v[i]); 
		
		u = fld->u[i];
		v = fld->v[i];
		sig = fld->sig[i];
		pres = fld->pres[i];
	
		
		omk = calc_omk(rc);
		
		fld->dtu[i] += I*m*omk*u + 2*(omk + omf) * v 
						- drpres/dbar[i] + sig*drpbar[i]/(dbar[i]*dbar[i]);
	
		fld->dtv[i] += I*m*omk*v - (2*omf + omk*(2+q))*u + (I*m/rc)*(pres/dbar[i]);
	
		fld->dtsig[i] += I*m*omk*sig - (1/rc) * (
						dbar[i] * u + u*rc*drdbar[i] + rc*dbar[i]*dru
						- I*m*(dbar[i] * v) );
	 
		fld->dtpres[i] +=  I*m*omk*pres - u*drpbar - (gam*pbar[i]/rc) * (
						 u + r dru - I*m*v );
	
	}
	


} 


double calc_omk(double r) {
	return om0*pow(r,q);
}