


void wavekillbc(Mode *fld,double dt)
{
	int i;
	double R,tau,x,dtdtau;
	const double x_inf = (fld->r[NG])*1.2;
	const double x_sup = (fld->Lr)*0.8;
	const double tau0 = .1
	
	for(i=istart;i<iend;i++) {
		x = fld->r[i]
		R=0;
		if (x > x_sup) R = (x-x_sup)/(fld->r[Nr+NG-1] - x_sup);
		if (x < x_inf) R = (x_inf - x)/(x_inf - fld->r[NG]);

		R *= R;
		tau = tau0;

		if (R>0.0) {
			tau /= R; 
			dtdtau = dt/tau;
			fld->u[i] = (fld->u[i])/(1+dtdtau );
			fld->v[i] = (fld->v[i])/(1+dtdtau );
			fld->sig[i] = (fld->sig[i])/(1+dtdtau);				

		}
	}
	return;
}