class Field():
	def __init__(self,t,outdir='',loadghost=False,loadrhs=False,rhsnum=0,NG=1,etol=1e-8):	
		outdir = check_dir(outdir)
		self.outdir = outdir
		disk = loadtxt(outdir+'disk.dat')
		self.r = disk[:,0]
		self.nlr = disk[:,1]
		self.nr = len(self.r)
		self.dr = diff(self.r)[0]
		self.width = self.nr*self.dr
		self.hor = disk[:,2]
		self.c2 = disk[:,3]
		self.nus = disk[:,4]
		self.nub = disk[:,5]
		self.nu = self.nus + self.nub
		self.Q = disk[:,6]
		dat = loadtxt(outdir + 'output_'+str(t)+'.dat')
		self.u = dat[:,2]+1j*dat[:,3]
		self.v = dat[:,4] + 1j*dat[:,5]
		self.sig = dat[:,6] + 1j*dat[:,7]
		self.vyb = dat[:,8]
		self.omk = dat[:,9]
		self.omk0 = pow(self.nlr,-1.5)
		self.dbar = dat[:,10]
		self.E = (2*self.v - 1j*self.u) / (2*self.vyb)
		self.kappa = sqrt(4*self.omk**2 + 2*self.omk * gradient(self.omk,self.dr))
		
		if not loadghost:
			self.r = self.r[NG:-NG]
			self.nlr = self.nlr[NG:-NG]
			self.hor = self.hor[NG:-NG]
			self.c2 = self.c2[NG:-NG]
			self.Q = self.Q[NG:-NG]
			self.nus = self.nus[NG:-NG]
			self.nub = self.nub[NG:-NG]
			self.nu = self.nu[NG:-NG]
			self.u = self.u[NG:-NG]
			self.v = self.v[NG:-NG]
			self.sig = self.sig[NG:-NG]
			self.vyb = self.vyb[NG:-NG]
			self.omk = self.omk[NG:-NG]
			self.omk0 = self.omk0[NG:-NG]
			self.dbar = self.dbar[NG:-NG]
			self.E = self.E[NG:-NG]
			self.kappa = self.kappa[NG:-NG]
			self.nr -= 2*NG
			
#		self.E = (abs(real(self.E))>etol).astype('int')*real(self.E) \
#				+1j*(abs(imag(self.E))>etol).astype('int') * imag(self.E)
#		self.E *= (((abs(real(self.E))<=etol).astype('int') * (abs(imag(self.E))<=etol).astype('int')+1)%2)
# 		for i in range(self.nr):
# 			if abs(real(self.E[i])) <= 1e-8:
# 				if abs(imag(self.E[i])) <= 1e-8:
# 					self.E[i] = 0
# 				else:
# 					self.E[i] = 1j*imag(self.E[i])
# 			else:
# 				if abs(imag(self.E[i])) <= 1e-8:
# 					self.E[i] = real(self.E[i])
		
		self.NG = NG
		if loadrhs:
			rhs = loadtxt('rhs_'+str(rhsnum)+'.dat')
			self.dtu = rhs[:,2]+rhs[:,3]*1j
			self.dtv = rhs[:,4] + rhs[:,5]*1j
			self.dts = rhs[:,6] + rhs[:,7]*1j
			if not loadghost:
				self.dtu = self.dtu[NG:-NG]
				self.dtv = self.dtv[NG:-NG]
				self.dts = self.dts[NG:-NG]
				
		else:
			self.dtu = 0
			self.dtv = 0
			self.dts = 0
		
		
		self.ind_dbar = polyfit(self.r,log(self.dbar),1)[0]
		self.ind_omk = polyfit(self.r,log(self.omk),1)[0]
		self.ind_c2 = polyfit(self.r,log(self.c2),1)[0]
		self.ind_hor = polyfit(self.r,log(self.hor),1)[0]
		self.ind_nu = polyfit(self.r,log(self.nu),1)[0]
		
		
		try:
			starpos = loadtxt(outdir + 'CentralStar.dat')
			try:
				(self.xstar,self.ystar,self.rstar,self.phistar) = starpos[t,1:]
			except:	
				(self.xstar,self.ystar,self.rstar,self.phistar)= (0,0,0,0)	
		except:
			(self.xstar,self.ystar,self.rstar,self.phistar)= (0,0,0,0)
			
		
		
	def plot(self,q,linestyle='-',logr=True,Nph=500,xnorm=0,ynorm=0,xlims=None,ylims=None,ncontours=20):
	
		if q not in ['all','u','v','sig','E','nu','c2','hor','Q','omk','dbar','vybar','dtu','dtv','dts','e','w','ex','ey'] \
		and q not in ['vr','vph','vphp','dens','densp']:
			print 'Not Valid Variable Name'
			return
		
		if q in ['vr','vph','vphp','dens','densp']:
			phi = linspace(-pi,pi,Nph)
			x = zeros((len(fld.nlr),Nph))
			y = zeros((len(fld.nlr),Nph))
			dat = zeros((len(fld.nlr),Nph))
			for i in range(Nph):
				x[:,i] = fld.nlr * cos(phi[i])
				y[:,i] = fld.nlr * sin(phi[i])
		
		
			
		if logr:
			r = copy(self.r)
			xname = '$\ln r$'
		else:
			r = copy(self.nlr)
			xname = '$r$'
			
		if xnorm != 0:
			r /= xnorm
			
		if q=='all':
			fig,((ax_ru,ax_rv,ax_rs,ax_re),(ax_iu,ax_iv,ax_is,ax_ie)) = subplots(2,4,sharex='col')
			
		
			ax_iu.set_xlabel(xname,fontsize='large')
			ax_iv.set_xlabel(xname,fontsize='large')
			ax_is.set_xlabel(xname,fontsize='large')
			ax_ie.set_xlabel(xname,fontsize='large')
		
			
			ax_ru.set_ylabel('$Re(u)$',fontsize='large')
			ax_rv.set_ylabel('$Re(v)$',fontsize='large')
			ax_rs.set_ylabel('$Re(\\sigma)$',fontsize='large')
			ax_re.set_ylabel('$e$',fontsize='large')
			
			
			ax_iu.set_ylabel('$Im(u)$',fontsize='large')
			ax_iv.set_ylabel('$Im(v)$',fontsize='large')
			ax_is.set_ylabel('$Im(\\sigma)$',fontsize='large')
			ax_ie.set_ylabel('$\\omega/\\pi$',fontsize='large')
			

			ax_ru.set_title('u',fontsize='large')
			ax_rv.set_title('v',fontsize='large')
			ax_rs.set_title('$\\sigma$',fontsize='large')
			ax_re.set_title('E',fontsize='large')


			ax_ru.plot(r,real(fld.u),label=r'$Re(u)$')
			ax_rv.plot(r,real(fld.v),label=r'$Re(v)$')
			ax_rs.plot(r,real(fld.sig),label=r'$Re(\\sigma)$')
			ax_re.plot(r,abs(fld.E),label=r'$e$')
		
			ax_iu.plot(r,imag(fld.u),label=r'$Im(u)$')
			ax_iv.plot(r,imag(fld.v),label=r'$Im(v)$')
			ax_is.plot(r,imag(fld.sig),label=r'$Im(\\sigma)$')
			ax_ie.plot(r,angle(fld.E)/pi,label=r'$\\omega/\\pi$')
	
		
		if q=='u':
			fig,(ax1,ax2)=subplots(2,sharex=True)
			ax1.set_title('u')
			ax1.plot(r,real(self.u),linestyle,label=r'$Re(u)$')
			ax1.set_ylabel('Re(u)')
			ax2.plot(r,imag(self.u),linestyle,label=r'$Im(u)$')
			ax2.set_ylabel('Im(u)')
			ax2.set_xlabel(xname)
		if q=='v':
			fig,(ax1,ax2)=subplots(2,sharex=True)
			ax1.set_title('v')
			ax1.plot(r,real(self.v),linestyle,label='$Re(v)$')
			ax1.set_ylabel('Re(v)')
			ax2.plot(r,imag(self.v),linestyle,label='$Im(v))$')
			ax2.set_ylabel('Im(v)')
			ax2.set_xlabel(xname)
		if q=='sig':
			fig,(ax1,ax2)=subplots(2,sharex=True)
			ax1.set_title('$\\sigma / <\\Sigma>$')
			ax1.plot(r,real(self.sig),linestyle,label='$Re(\\sigma))$')
			ax1.set_ylabel('Re($\\sigma$)')
			ax2.plot(r,imag(self.sig),linestyle,label='$Im(\\sigma))$')
			ax2.set_ylabel('Im($\\sigma$)')
			ax2.set_xlabel(xname)
		if q=='vybar':
			fig,ax = subplots()
			ax.set_title('$<v_\\phi>$')
			ax.plot(r,self.vyb,linestyle)
			ax.set_xlabel(xname)
			
		if q=='dbar':
			fig,ax = subplots()
			ax.set_title('$<\\Sigma>$')
			ax.plot(r,self.dbar,linestyle)
			ax.set_xlabel(xname)
			
		if q=='omk':
			fig,ax = subplots()
			ax.set_title('$\\Omega_k$')
			ax.plot(r,self.omk,linestyle)
			ax.set_xlabel(xname)
		if q=='nu':
			fig,(ax1,ax2,ax3)= subplots(3,sharex=True)
			ax1.set_title('$\\nu$')
			ax1.set_ylabel('$\\nu$')
			ax2.set_ylabel('$\\nu_s$')
			ax3.set_ylabel('$\\nu_b$')
			ax1.plot(r,self.nu,linestyle)
			ax2.plot(r,self.nus,linestyle)
			ax3.plot(r,self.nub,linestyle)
			ax3.set_xlabel(xname)
		if q=='c2':
			fig,ax = subplots()
			ax.set_title('$c^2$')
			ax.plot(r,self.c2,linestyle)
			ax.set_xlabel(xname)
		if q=='hor':	
			fig,ax = subplots()
			ax.set_title('$h/r$')
			ax.plot(r,self.hor,linestyle)
			ax.set_xlabel(xname)
			
		if q=='Q':	
			fig,ax = subplots()
			ax.set_title('$Q$')
			ax.plot(r,self.Q,linestyle)
			ax.set_xlabel(xname)
			
		if q=='E':
			fig,(ax1,ax2,ax3,ax4)=subplots(4,sharex=True)
			ax1.set_title('Eccentricity')
			ax1.plot(r,real(self.E),linestyle)
			ax1.set_ylabel('$e_x$')
			ax2.plot(r,imag(self.E),linestyle)
			ax2.set_ylabel('$e_y$')
			ax3.plot(r,abs(self.E),linestyle)
			ax3.set_ylabel('e')
			ax4.plot(r,angle(self.E)/pi,linestyle)
			ax4.set_ylabel('$\\omega/\\pi$')
			ax4.set_xlabel(xname)
		
		if q=='e':
			fig,ax = subplots()
			ax.set_title('$e$')
			if ynorm == 0:
				ax.plot(r,real(abs(self.E)),linestyle)
			else:
				ax.plot(r,real(abs(self.E))/ynorm,linestyle)
			ax.set_xlabel(xname)
		if q=='w':
			fig,ax = subplots()
			ax.set_title('$\\omega/\\pi$')
			ax.plot(r,angle(self.E)/pi,linestyle)
			ax.set_xlabel(xname)
		if q=='ex':
			fig,ax = subplots()
			ax.set_title('$e_x$')
			ax.plot(r,real(self.E),linestyle)
			ax.set_xlabel(xname)	
		if q=='ey':
			fig,ax = subplots()
			ax.set_title('$e_y$')
			ax.plot(r,imag(self.E),linestyle)
			ax.set_xlabel(xname)	
		if q=='dtu':
			fig,(ax1,ax2)=subplots(2,sharex=True)
			ax1.set_title('dtu')
			ax1.plot(r[self.NG:-self.NG],real(self.dtu),linestyle,label=r'$Re(dtu)$')
			ax1.set_ylabel('Re(dtu)')
			ax2.plot(r[self.NG:-self.NG],imag(self.dtu),linestyle,label=r'$Im(dtu)$')
			ax2.set_ylabel('Im(dtu)')
			ax2.set_xlabel(xname)
		if q=='dtv':
			fig,(ax1,ax2)=subplots(2,sharex=True)
			ax1.set_title('dtv')
			ax1.plot(r[self.NG:-self.NG],real(self.dtv),linestyle,label='$Re(dtv)$')
			ax1.set_ylabel('Re(dtv)')
			ax2.plot(r[self.NG:-self.NG],imag(self.dtv),linestyle,label='$Im(dtv))$')
			ax2.set_ylabel('Im(dtv)')
			ax2.set_xlabel(xname)
		if q=='dts':
			fig,(ax1,ax2)=subplots(2,sharex=True)
			ax1.set_title('dts')
			ax1.plot(r[self.NG:-self.NG],real(self.dts),linestyle,label='$Re(dts)$')
			ax1.set_ylabel('Re(dts)')
			ax2.plot(r[self.NG:-self.NG],imag(self.dts),linestyle,label='$Im(dts)$')
			ax2.set_ylabel('Im(dts)')	
			ax2.set_xlabel(xname)
		
		if q in ['vr','vph','vphp','dens','densp']:
			if q in ['dens','densp']:	
				for i in range(Nph):
					dat[:,i] = self.dbar*(1 + 2*real(self.sig*exp(-1j*phi[i])))
					if q=='dens':
						tstr='$<\\Sigma>$'
					else:
						dat[:,i] = dat[:,i]/self.dbar - 1;	
						tstr='$\\Delta \\Sigma / <	\\Sigma > $'
						
			if q == 'vr':
				for i in range(Nph):
					dat[:,i] = 2*real(self.u*exp(-1j*phi[i]))
				tstr='$v_r$'
			if q == 'vph':
				for i in range(Nph):
					dat[:,i] = 2*real(self.v*exp(-1j*phi[i])) + self.vyb
				tstr='$v_\\phi$'
			if q == 'vphp':
				for i in range(Nph):
					dat[:,i] = 2*real(self.v*exp(-1j*phi[i]))
				tstr='$\\Delta v_\\phi$'

			figure()
			pcolormesh(x,y,dat)
			colorbar()
			contour(x,y,dat,ncontours,colors='k')
			plot(self.xstar,self.ystar,'k*',markersize=10)
			plot(0,0,'k.')
			
			if ylims != None:
				ylim(ylims)
			title(tstr)
			xlabel('x')
			ylabel('y')
		
	
		if xlims != None:
				xlim(xlims)
		show()
		
		
		
	def disk_summary(self,logr=True):
		fig,((ax_omk,ax_c2,ax_nu),(ax_sig,ax_hor,ax_Q)) = subplots(2,3,sharex='col')
			
		if logr:
			xstr =  '$\ln r$'
			r = fld.r
		else:
			xstr = '$r$'
			r = fld.nlr
		
		ax_sig.set_xlabel(xstr,fontsize='large')
		ax_hor.set_xlabel(xstr,fontsize='large')
		ax_Q.set_xlabel(xstr,fontsize='large')
		
		ax_omk.set_title('$\\Omega$',fontsize='large')
		ax_c2.set_title('$c^2$',fontsize='large')
		ax_nu.set_title('$\\log_{10}(\\nu)$',fontsize='large')
		
		ax_sig.set_title('$<\\Sigma>$',fontsize='large')
		ax_hor.set_title('$H/r$',fontsize='large')
		ax_Q.set_title('Toomre $Q$',fontsize='large')
		
		
		omk_str = '$d\ln \\Omega / d\ln r$ = %.1f' % fld.ind_omk
		sig_str = '$d\ln <\\Sigma> / d\ln r$ = %.1f' % fld.ind_dbar
		c2_str = '$d\ln c^2 / d\ln r$ = %.1f' % fld.ind_c2
		hor_str = '$d\ln ( H/r ) / d\ln r$ = %.1f' % fld.ind_hor
		nu_str = '$d\ln \\nu / d\ln r$ = %.1f' % fld.ind_nu
		
		ax_omk.plot(r,fld.omk,label=omk_str)
		ax_c2.plot(r,fld.c2,label=c2_str)
		ax_nu.plot(r,log10(fld.nu),label=nu_str)
		
		ax_sig.plot(r,fld.dbar,label=sig_str)
		ax_hor.plot(r,fld.hor,label=hor_str)
		ax_Q.plot(r,fld.Q)
		
		ax_omk.legend(loc='best')
		ax_c2.legend(loc='best')
		ax_nu.legend(loc='best')
		ax_sig.legend(loc='best')
		ax_hor.legend(loc='best')
		
		
		show()
		
	def plot_cfl(self,logr=True):
		
		if logr:
			r = fld.r
			xstr = '$\ln r$'
		else:
			r = fld.nlr
			xstr = '$r$'
			
			
		dtc = fld.dr * fld.nlr / sqrt(fld.c2)
		
		dtnu = (fld.dr * fld.nlr)**2 / fld.nu
		
		
		figure()
		xlabel(xstr,fontsize='large')
		ylabel('$\\log_{10}(\\Delta t)$',fontsize='large')
		plot(r,log10(dtc),label=r'Sound Crossing Time, $t_{total}$ = %.2e' % sum(dtc))
		plot(r,log10(dtnu),label=r'Viscous Crossing Time, $t_{total}$ = %.2e' % sum(dtnu))
		legend(loc='best')
		return
	
			
	def draw_ellipse(self,num_ellipse,(xc,yc)=(0,0),Nph=500):
		pgrid = linspace(0,2*pi,Nph)
		# (Np , Nr)		
		
		
		x = zeros((Nph,self.nr))
		y = zeros((Nph,self.nr))
		for i in range(Nph):
			# semi latus rectum
			l = self.nlr * (self.vyb  + real(self.v*exp(-1j*pgrid[i])))
			l0 = self.nlr*self.nlr * self.omk0
			p = self.nlr*(l/l0)*(l/l0)
			theta = pgrid[i] - angle(self.E)
			a = p/(1-abs(self.E)**2)
			b = p/sqrt(1-abs(self.E)**2)
			x[i,:] = xc +a*cos(theta) 
			y[i,:] = yc + b*sin(theta)
			
			
		figure();	
		size_x = abs(x).max()
		size_y = abs(y).max()
		xlim((-size_x,size_x))
		ylim((-size_y,size_y))
		plot(self.xstar,self.ystar,'k*',markersize=10)
		for i in range(self.nr)[::self.nr/num_ellipse]:
			plot(x[:,i],y[:,i],'-k')
				
		return x,y
	
	def plotdens(self,Nph=500,starpos=(0,0)):
		phi = linspace(-pi,pi,Nph)
		dens = zeros((len(fld.nlr),Nph))	
		densp = zeros((len(fld.nlr),Nph))
		x = zeros((len(fld.nlr),Nph))
		y = zeros((len(fld.nlr),Nph))
		for i in range(Nph):
			dens[:,i] = fld.dbar*(1 + 2*real(fld.sig*exp(-1j*phi[i])));
			densp[:,i] = dens[:,i]/fld.dbar - 1;
			x[:,i] = fld.nlr * cos(phi[i])
			y[:,i] = fld.nlr * sin(phi[i])
		
		figure()
		pcolormesh(x,y,densp); colorbar()
		plot([starpos[0]],[starpos[1]],'*')
		title('$\\frac{\Sigma - <\Sigma>}{<\Sigma>}$')
		return
	

			
				
	
def animate(q,t,dt=1,linestyle='-',dat=None,fld0=None,logr=True):
	if q not in ['u','v','sig','E','e','w','ex','ey']:
		print 'Not Valid Variable Name'
		return
		
	if fld0==None:
		fld0 = Field(0)
		
	if dat==None:
		dat = zeros((fld0.nr,len(t)),dtype='complex')
		for i,j in enumerate(t):
			print 'Loading time ',j
			fld = Field(j)
			if q=='u':
				dat[:,i] = fld.u
			if q=='v':
				dat[:,i] = fld.v
			if q=='sig':
				dat[:,i] = fld.sig
			if q=='E':
				dat[:,i] = fld.E
			if q=='e':
				dat[:,i] = real(abs(fld.E))
			if q=='w':
				dat[:,i] = angle(fld.E)
			if q=='ex':
				dat[:,i] = real(fld.E)
			if q=='ey':
				dat[:,i] = imag(fld.E)				
	
	if logr:
		r = fld0.r
		xstr = '$\ln r$'
	else:
		r = fld0.nlr
		xstr = 'r'
	
	if q in ['e','w','ex','ey']:
		fig = figure()
		ax = fig.add_subplot(111)
		ax.set_title(q)
		ax.set_xlabel(xstr)
		l, = ax.plot(r,dat[:,0])
		lrange = (dat.min(),dat.max())
		for i in range(dat.shape[1]):
			ax.set_title(q + ', t='+str(t[i]*dt))
			l.set_ydata(dat[:,i])
			ax.set_ylim(lrange)
			fig.canvas.draw()
	
		
	if q=='u':
		fig = figure()
		ax1=fig.add_subplot(211)
		ax2 = fig.add_subplot(212,sharex=ax1)
		ax1.set_title('u')
		ax1.set_ylabel('Re(u)')
		ax2.set_ylabel('Im(u)')
		ax2.set_xlabel(xstr)
		l1,= ax1.plot(r,real(fld0.u),linestyle)
		l2,= ax2.plot(r,imag(fld0.u),linestyle)
		l1range = (real(dat).min(),real(dat).max())
		l2range = (imag(dat).min(),imag(dat).max())
		for i in range(dat.shape[1]):
			ax1.set_title('u, t='+str(t[i]*dt))
			l1.set_ydata(real(dat[:,i]))
			l2.set_ydata(imag(dat[:,i]))
			ax1.set_ylim(l1range)
			ax2.set_ylim(l2range)
			fig.canvas.draw()
		
			
		
	if q=='v':
		fig = figure()
		ax1=fig.add_subplot(211)
		ax2 = fig.add_subplot(212,sharex=ax1)
		ax1.set_title('v, t=0')
		ax1.set_ylabel('Re(v)')
		ax2.set_ylabel('Im(v)')
		ax2.set_xlabel(xstr)
		l1,= ax1.plot(r,real(fld0.v),linestyle)
		l2,= ax2.plot(r,imag(fld0.v),linestyle)
		l1range = (real(dat).min(),real(dat).max())
		l2range = (imag(dat).min(),imag(dat).max())
		for i in range(dat.shape[1]):
			ax1.set_title('v, t='+str(t[i]*dt))
			l1.set_ydata(real(dat[:,i]))
			l2.set_ydata(imag(dat[:,i]))
			ax1.set_ylim(l1range)
			ax2.set_ylim(l2range)
			fig.canvas.draw()
	if q=='sig':
		fig = figure()
		ax1=fig.add_subplot(211)
		ax2 = fig.add_subplot(212,sharex=ax1)
		
		ax1.set_title('$\\sigma / <\\Sigma>$, t=0')
		ax1.set_ylabel('Re($\\sigma$)')
		ax2.set_ylabel('Im($\\sigma$)')
		ax2.set_xlabel(xstr)
		l1,= ax1.plot(r,real(fld0.sig),linestyle)
		l2,= ax2.plot(r,imag(fld0.sig),linestyle)
		
		l1range = (real(dat).min(),real(dat).max())
		l2range = (imag(dat).min(),imag(dat).max())
		fig.canvas.draw()
		for i in range(dat.shape[1]):
			ax1.set_title('$\\sigma / <\\Sigma>$, t='+str(t[i]*dt))
			l1.set_ydata(real(dat[:,i]))
			l2.set_ydata(imag(dat[:,i]))
			ax1.set_ylim(l1range)
			ax2.set_ylim(l2range)
			fig.canvas.draw()
	if q=='E':
		fig = figure()
		ax1=fig.add_subplot(411)
		ax2 = fig.add_subplot(412,sharex=ax1)
		ax3 = fig.add_subplot(413,sharex=ax1)
		ax4 = fig.add_subplot(414,sharex=ax1)
		ax1.set_title('Eccentricity, t=0')
		ax1.set_ylabel('$e_x$')
		ax2.set_ylabel('$e_y$')
		ax3.set_ylabel('e')
		ax4.set_ylabel('$\\omega/\\pi$')
		ax4.set_xlabel(xstr)	
		
		l1,= ax1.plot(r,real(fld0.E),linestyle)
		l2,= ax2.plot(r,imag(fld0.E),linestyle)
		l3,= ax3.plot(r,abs(fld0.E),linestyle)
		l4,= ax4.plot(r,angle(fld0.E)/pi,linestyle)
		
		l1range = (real(dat).min(),real(dat).max())
		l2range = (imag(dat).min(),imag(dat).max())
		l3range = (abs(dat).min(),abs(dat).max())
		l4range = (-1,1)
		fig.canvas.draw()
		for i in range(len(t)):
			ax1.set_title('Eccentricity, t='+str(t[i]*dt))
			l1.set_ydata(real(dat[:,i]))
			l2.set_ydata(imag(dat[:,i]))
			l3.set_ydata(abs(dat[:,i]))
			l4.set_ydata(angle(dat[:,i])/pi)

			ax1.set_ylim(l1range)
			ax2.set_ylim(l2range)
			ax3.set_ylim(l3range)
			ax4.set_ylim(l4range)			
			
			fig.canvas.draw()
	
	
					
	return dat,fld0		
	

def animate_ellipse(t,num_ellipse,(xc,yc)=(0,0),Nph=500):
		
		
	if fld0==None:
		fld0 = Field(0)
	
	
	if dat==None:
		pgrid=linspace(0,2*pi,Nph)
		ind = fld0.r==fld0.r[::fld.nr/num_ellipse]
		l0 = fld0.nlr[ind]*fld0.nlr[ind] * fld0.omk0[ind]
		x = zeros((Nph,len(l0),len(t)))
		y = zeros((Nph,len(l0),len(t)))
		
		for i,j in enumerate(t):
			print 'Loading time', j
			fld=Field(j)
			w = angle(fld.E[ind])
			
			l = fld.nlr[ind] * (fld.vyb[ind]  + real(fld.v[ind]*exp(-1j*pgrid[i])))
			p = fld.nlr[ind]*(l/l0)*(l/l0)
			theta = pgrid[i] - angle(fld.E[ind])
			a = p/(1-abs(fld.E)**2)
			b = p/sqrt(1-abs(self.E)**2)
			x[i,:] = xc +a*cos(theta) 
			y[i,:] = yc + b*sin(theta)			
	
	if q=='u':
		fig = figure()
		ax1=fig.add_subplot(211)
		ax2 = fig.add_subplot(212,sharex=ax1)
		ax1.set_title('u')
		ax1.set_ylabel('Re(u)')
		ax2.set_ylabel('Im(u)')
		ax2.set_xlabel('$\ln r$')
		l1,= ax1.plot(fld0.r,real(fld0.u),linestyle)
		l2,= ax2.plot(fld0.r,imag(fld0.u),linestyle)
		l1range = (real(dat).min(),real(dat).max())
		l2range = (imag(dat).min(),imag(dat).max())
		for i in range(dat.shape[1]):
			ax1.set_title('v, t='+str(t[i]*dt))
			l1.set_ydata(real(dat[:,i]))
			l2.set_ydata(imag(dat[:,i]))
			ax1.set_ylim(l1range)
			ax2.set_ylim(l2range)
			fig.canvas.draw()


def plotstar(outdir='',inner_rad=None,linestyle='-*'):
		
	dat = loadtxt(outdir + 'CentralStar.dat')
	
	growth_rate = polyfit(dat[1:,0]/(2*pi),log(dat[1:,3]),1)
	print 'Growth rate is ', growth_rate
	glabel = '$d\ln r_\star / dt$ = %.2e' % growth_rate[0]
	
	fig = figure()
	ax = fig.add_subplot(121)
	ax.set_title('Star Position',fontsize='large')
	ax.set_xlabel('x',fontsize='large')
	ax.set_ylabel('y',fontsize='large')
	ax.plot(dat[:,1],dat[:,2],linestyle)
	if inner_rad != None:
		if inner_rad <= 10*dat[:,3].max():
			ax.plot(inner_rad*cos(linspace(0,2*pi,100)), inner_rad*sin(linspace(0,2*pi,100)),'-r')
		else:
			print 'Disk inner radius too far out!'
	ax2 = fig.add_subplot(122)
	ax2.set_title('Star radius',fontsize='large')
	ax2.set_xlabel('$ t_{orb}$',fontsize='large')
	ax2.set_ylabel('$ \ln r_\star $',fontsize='large')
	ax2.plot(dat[:,0]/(2*pi),log(dat[:,3]),linestyle,label=glabel)
	if inner_rad != None:
		if inner_rad <= 10*dat[:,3].max():
			ax2.plot(dat[:,0]/(2*pi),log(inner_rad)*ones(len(dat[:,0])),'-r')
	ax2.legend(loc='lower right')
	
def calc_growth_rate(times, dt, outputdir=''):
	outputdir = check_dir(outputdir)
 
	ecc = zeros((len(times),1))
	print 'Loading times'
	for i,t in enumerate(times):
		fld=Field(t,outdir=outputdir)
		ecc[i] = real(abs(fld.E)).max(axis=0)
	ecc = log(ecc / ecc.max())
	
	t = times * dt / (2*pi)
	
	growth = polyfit(t,ecc,1)
	estr = '$d\ln e / dt$ = %.2e' % growth[0][0]
	
	fig = figure()
	plot(t,ecc)
	plot(t, growth[0]*t + growth[1],'--',label=estr)
	xlabel('t',fontsize='large')
	ylabel('$\ln ( e_{peak}/ e_{max} )$',fontsize='large')
	legend(loc='best')
	return growth, ecc, t
def check_dir(directory):
	if len(directory) != 0:
		if directory[-1] != '/':
			directory += '/'
	return directory
			
def make_star_eccen_plots(none_dir,ind_dir,sg_dir,both_dir,tend):
	sg_dir = check_dir(sg_dir)
	ind_dir = check_dir(ind_dir)
	both_dir = check_dir(both_dir)
	none_dir = check_dir(none_dir)
	
	dat1 = loadtxt(ind_dir + 'CentralStar.dat')
	dat3 = loadtxt(both_dir + 'CentralStar.dat')
	
	
	dat1 = dat1[:tend,:]
	dat3 = dat3[:tend,:]
	
	t = dat1[:,0]
	
	fld0 = Field(0,outdir=none_dir)
	fld1 = Field(0,outdir=ind_dir)
	fld2 = Field(0,outdir=sg_dir)
	fld3 = Field(0,outdir=both_dir)
	
	t /= (2*pi/fld0.omk[0])
	
	ecc0 = zeros((len(t),1))
	ecc1 = zeros((len(t),1))
	ecc2 = zeros((len(t),1))
	ecc3 = zeros((len(t),1))
	
	print 'Loading times'
	for i in range(tend):
	
		fld=Field(i,outdir=none_dir)
		ecc0[i] = real(abs(fld.E)).max(axis=0)
		
		fld=Field(i,outdir=ind_dir)
		ecc1[i] = real(abs(fld.E)).max(axis=0)
		
		fld=Field(i,outdir=sg_dir)
		ecc2[i] = real(abs(fld.E)).max(axis=0)
		
		fld=Field(i,outdir=both_dir)
		ecc3[i] = real(abs(fld.E)).max(axis=0)
	print 'Done'	
	
	ecc0 = log(ecc0/ecc0.max())
	ecc1 = log(ecc1/ecc1.max())
	ecc2 = log(ecc2/ecc2.max())
	ecc3 = log(ecc3/ecc3.max())
	
	growth0 = polyfit(t,ecc0,1) 
	growth1 = polyfit(t,ecc1,1) 
	growth2 = polyfit(t,ecc2,1) 
	growth3 = polyfit(t,ecc3,1) 
	
	estr0 = '$d\ln e / dt$ = %.2e' % growth0[0][0]
	estr1 = '$d\ln e / dt$ = %.2e' % growth1[0][0]
	estr2 = '$d\ln e / dt$ = %.2e' % growth2[0][0]
	estr3 = '$d\ln e / dt$ = %.2e' % growth3[0][0]
	
	
	x0,y0,r0 = (0,0,0)
	x2,y2,r2 = (0,0,0)
	x1 = dat1[:,1]
	y1 = dat1[:,2]
	r1 = log(dat1[:,3] / fld1.nlr[0])
	x3 = dat3[:,1]
	y3 = dat3[:,2]
	r3 = log(dat3[:,3] / fld3.nlr[0])
	
	
	rgrowth1 = polyfit(t,ecc1,1) 
	rgrowth3 = polyfit(t,ecc3,1) 
	
	rstr1 = '$d\ln r_\star / dt$ = %.2e' % rgrowth1[0][0]
	rstr3 = '$d\ln r_\star / dt$ = %.2e' % rgrowth3[0][0]
	
	
	fig = figure()
	ax_none_star = fig.add_subplot(431)
	ax_none_r = fig.add_subplot(432)
	ax_none_e = fig.add_subplot(433)
	
	ax_ind_star = fig.add_subplot(434)
	ax_ind_r = fig.add_subplot(435)
	ax_ind_e = fig.add_subplot(436)
	
	ax_sg_star = fig.add_subplot(437)
	ax_sg_r = fig.add_subplot(438)
	ax_sg_e = fig.add_subplot(439)
	
	ax_both_star = fig.add_subplot(4,3,10)
	ax_both_r = fig.add_subplot(4,3,11)
	ax_both_e = fig.add_subplot(4,3,12)
	
	
	
	ax_none_star.set_title('Star Position',fontsize='large')
	ax_none_r.set_title('Star Radius over Time',fontsize='large')
	ax_none_e.set_title('Maximum Disk Eccentricity over Time',fontsize='large')
	
	ax_both_star.set_xlabel('$x$',fontsize='large')
	ax_both_r.set_xlabel('$t/P_{inner}$',fontsize='large')
	ax_both_e.set_xlabel('$t/P_{inner}$',fontsize='large')
	
	ax_none_star.set_ylabel('$y$',fontsize='large')
	ax_ind_star.set_ylabel('$y$',fontsize='large')
	ax_sg_star.set_ylabel('$y$',fontsize='large')
	ax_both_star.set_ylabel('$y$',fontsize='large')
	
	
	ax_none_r.set_ylabel('$\ln ( r_\star / r_i )$',fontsize='large')
	ax_ind_r.set_ylabel('$\ln ( r_\star / r_i )$',fontsize='large')
	ax_sg_r.set_ylabel('$\ln ( r_\star / r_i )$',fontsize='large')
	ax_both_r.set_ylabel('$\ln ( r_\star / r_i )$',fontsize='large')
	
	ax_none_e.set_ylabel('$\ln ( e_{peak}/ e_{max} )$',fontsize='large')
	ax_ind_e.set_ylabel('$\ln ( e_{peak}/ e_{max} )$',fontsize='large')
	ax_sg_e.set_ylabel('$\ln ( e_{peak}/ e_{max} )$',fontsize='large')
	ax_both_e.set_ylabel('$\ln ( e_{peak}/ e_{max} )$',fontsize='large')
	
	ax_none_star.plot(x0,y0,'-*')
	ax_ind_star.plot(x1,y1,'-*')
	ax_sg_star.plot(x2,y2,'-*')
	ax_both_star.plot(x3,y3,'-*')
	
	ax_ind_r.plot(t,r1,'-*',label=rstr1)
	ax_both_r.plot(t,r3,'-*',label=rstr3)
	
	ax_none_e.plot(t,ecc0,label=estr0)
	ax_ind_e.plot(t,ecc1,label=estr1)
	ax_sg_e.plot(t,ecc2,label=estr2)
	ax_both_e.plot(t,ecc3,label=estr3)
	
	
	
	ind_star_max = max( ax_ind_star.get_xlim() + ax_ind_star.get_ylim())
	ind_star_min = min( ax_ind_star.get_xlim() + ax_ind_star.get_ylim())
	ind_star_lim = max((ind_star_max,abs(ind_star_min)))
	ax_ind_star.set_xlim( (-ind_star_lim, ind_star_lim) )
	ax_ind_star.set_ylim( (-ind_star_lim, ind_star_lim) )
	
	both_star_max = max( ax_both_star.get_xlim() + ax_both_star.get_ylim())
	both_star_min = min( ax_both_star.get_xlim() + ax_both_star.get_ylim())
	both_star_lim = max((both_star_max,abs(both_star_min)))
	ax_both_star.set_xlim( (-both_star_lim, both_star_lim) )
	ax_both_star.set_ylim( (-both_star_lim, both_star_lim) )
	
	
	ax_none_e.legend(loc='upper right')
	
	ax_ind_r.legend(loc='upper right')
	ax_ind_e.legend(loc='upper right')
	
	ax_sg_e.legend(loc='lower right')
	
	ax_both_r.legend(loc='lower right')
	ax_both_e.legend(loc='lower right')
	
	mean(ax_none_e.get_xlim())
	
	x_coord_label_none = ax_none_e.get_xlim()[1]
	x_coord_label_ind = ax_ind_e.get_xlim()[1]
	x_coord_label_sg = ax_sg_e.get_xlim()[1]
	x_coord_label_both = ax_both_e.get_xlim()[1]
	 
	y_coord_label_none = mean(ax_none_e.get_ylim()) 
	y_coord_label_ind = mean(ax_ind_e.get_ylim()) 
	y_coord_label_sg = mean(ax_sg_e.get_ylim()) 
	y_coord_label_both = mean(ax_both_e.get_ylim()) 
	
	bbox_props = dict(boxstyle='square',facecolor='white')
	ax_none_e.text(x_coord_label_none,y_coord_label_none,'Neither',fontsize=15,bbox=bbox_props)
	ax_ind_e.text(x_coord_label_ind,y_coord_label_ind,'Indirect only',fontsize=15,bbox=bbox_props)
	ax_sg_e.text(x_coord_label_sg,y_coord_label_sg,'SG only',fontsize=15,bbox=bbox_props)
	ax_both_e.text(x_coord_label_both,y_coord_label_both,'Both',fontsize=15,bbox=bbox_props)

	
	show()
	
		
def animate_real(q,t,xlims=None,ylims=None,Nph=500,ncontours=20):

	if q not in ['dens','densp','vr','vph','vphp','E']:
		print 'Not Valid Variable Name'
		return 

	fld0=Field(t[0])
	phi = linspace(-pi,pi,Nph)
	dat = zeros((len(fld0.nlr),Nph,len(t)))	

	x = zeros((len(fld0.nlr),Nph))
	y = zeros((len(fld0.nlr),Nph))
	
	for p in range(Nph):
		x[:,p] = fld0.nlr * cos(phi[p])
		y[:,p] = fld0.nlr * sin(phi[p])
	
	
	if q=='vr':
		tstr = '$v_r$'
	if q=='vphp':
		tstr = '$v_\phi - < v_\phi >$'
	if q=='vph':
		tstr = '$v_\phi$'
	if q=='dens':
		tstr = '$\\Sigma$'
	if q=='densp':
		tstr = '$(\\Sigma- < \\Sigma >)/<\\Sigma>$'
	if q=='E':
		tstr = '$e$'
	
	for i,j in enumerate(t):
		print 'Loading t = ' + str(i)
		fld=Field(j)
		for p in range(Nph):
			if q=='vr':
				dat[:,p,i] =  2*real(fld.u*exp(-1j*phi[p]))
			if q=='vphp':
				dat[:,p,i] =   2*real(fld.v*exp(-1j*phi[p]))
			if q=='vph':
				dat[:,p,i] =   2*real(fld.v*exp(-1j*phi[p])) + fld.vyb
			if q=='dens':
				dat[:,p,i] =  fld.dbar*(1 + 2*real(fld.sig*exp(-1j*phi[p])))
			if q=='densp':
				dat[:,p,i] =  fld.dbar*(1 + 2*real(fld.sig*exp(-1j*phi[p])))/fld.dbar - 1
			if q=='E':
				dat[:,p,i] =  2*abs(real(fld.E*exp(-1j*phi[p])))
	
	fig=figure()
	for i,j in enumerate(t[1:]):
		fig.clear()
		pcolormesh(x,y,dat[:,:,i])
		colorbar()
		contour(x,y,dat[:,:,i],ncontours,colors='k')
		title(tstr + ',    t = ' + str(j))
		if xlims != None:
			xlim(xlims)
		if ylims != None:
			ylim(ylims)
		fig.canvas.draw()
		
	
	
	return dat,x,y,fld0


def compare(q,fld_list,logr=True,linestyle='-'):

	if q not in ['all','u','v','sig','E','nu','c2','hor','omk','dbar','vybar','e','w','ex','ey']:
		print 'Not Valid Variable Name'
		return
		
		
	if logr:
		r = [fld.r for fld in fld_list]
		xname = '$\ln r$'
	else:
		r = [fld.nlr for fld in fld_list]		
		xname = '$r$'
	
	
	leg_str = [ str(i) for i in range(len(fld_list))]
	
	
	if q=='all':
	
		fig,((ax_ru,ax_rv,ax_rs,ax_re),(ax_iu,ax_iv,ax_is,ax_ie)) = subplots(2,4,sharex='col')
			
		
		ax_iu.set_xlabel(xname,fontsize='large')
		ax_iv.set_xlabel(xname,fontsize='large')
		ax_is.set_xlabel(xname,fontsize='large')
		ax_ie.set_xlabel(xname,fontsize='large')
	
		
		ax_ru.set_ylabel('$Re(u)$',fontsize='large')
		ax_rv.set_ylabel('$Re(v)$',fontsize='large')
		ax_rs.set_ylabel('$Re(\\sigma)$',fontsize='large')
		ax_re.set_ylabel('$e$',fontsize='large')
		
		
		ax_iu.set_ylabel('$Im(u)$',fontsize='large')
		ax_iv.set_ylabel('$Im(v)$',fontsize='large')
		ax_is.set_ylabel('$Im(\\sigma)$',fontsize='large')
		ax_ie.set_ylabel('$\\omega/\\pi$',fontsize='large')
		

		ax_ru.set_title('u',fontsize='large')
		ax_rv.set_title('v',fontsize='large')
		ax_rs.set_title('$\\sigma$',fontsize='large')
		ax_re.set_title('E',fontsize='large')



		for i,fld in enumerate(fld_list):
			ax_ru.plot(r[i],real(fld.u),label=r'$Re(u)$')
			ax_rv.plot(r[i],real(fld.v),label=r'$Re(v)$')
			ax_rs.plot(r[i],real(fld.sig),label=r'$Re(\\sigma)$')
			ax_re.plot(r[i],abs(fld.E),label=r'$e$')
	
			ax_iu.plot(r[i],imag(fld.u),label=r'$Im(u)$')
			ax_iv.plot(r,imag(fld.v),label=r'$Im(v)$')
			ax_is.plot(r,imag(fld.sig),label=r'$Im(\\sigma)$')
			ax_ie.plot(r,angle(fld.E)/pi,label=r'$\\omega/\\pi$')	
		
			
		
	if q=='u':
		fig,(ax1,ax2)=subplots(2,sharex=True)
		ax1.set_title('u')
		ax1.set_ylabel('Re(u)')
		ax2.set_ylabel('Im(u)')
		ax2.set_xlabel(xname)
		for i,fld in enumerate(fld_list):
			ax1.plot(r[i],real(fld.u),linestyle)
			ax2.plot(r[i],imag(fld.u),linestyle)
		
		
	if q=='v':
		fig,(ax1,ax2)=subplots(2,sharex=True)
		ax1.set_title('v')
		ax1.set_ylabel('Re(v)')
		ax2.set_ylabel('Im(v)')
		ax2.set_xlabel(xname)
		for i,fld in enumerate(fld_list):
			ax1.plot(r[i],real(fld.v),linestyle)
			ax2.plot(r[i],imag(fld.v),linestyle)
			
	if q=='sig':
		fig,(ax1,ax2)=subplots(2,sharex=True)
		ax1.set_title('$\\sigma / <\\Sigma>$')
		ax1.set_ylabel('Re($\\sigma$)')
		ax2.set_ylabel('Im($\\sigma$)')
		ax2.set_xlabel(xname)
		for i,fld in enumerate(fld_list):
			ax1.plot(r[i],real(fld.sig),linestyle)
			ax2.plot(r[i],imag(fld.sig),linestyle)
	if q=='vybar':
		fig,ax = subplots()
		ax.set_title('$<v_\\phi>$')
		ax.set_xlabel(xname)
		for i,fld in enumerate(fld_list):
			ax.plot(r[i],fld.vyb,linestyle)
		
	if q=='dbar':
		fig,ax = subplots()
		ax.set_title('$<\\Sigma>$')
		ax.set_xlabel(xname)
		for i,fld in enumerate(fld_list):
			ax.plot(r[i],fld.dbar,linestyle)
		
	if q=='omk':
		fig,ax = subplots()
		ax.set_title('$\\Omega_k$')
		ax.set_xlabel(xname)
		for i,fld in enumerate(fld_list):
			ax.plot(r[i],fld.omk,linestyle)
			
	if q=='nu':
		fig,(ax1,ax2,ax3)= subplots(3,sharex=True)
		ax1.set_title('$\\nu$')
		ax1.set_ylabel('$\\nu$')
		ax2.set_ylabel('$\\nu_s$')
		ax3.set_ylabel('$\\nu_b$')
		ax3.set_xlabel(xname)
		for i,fld in enumerate(fld_list):
			ax1.plot(r[i],fld.nu,linestyle)
			ax2.plot(r[i],fld.nus,linestyle)
			ax3.plot(r[i],fld.nub,linestyle)
			
	if q=='c2':
		fig,ax = subplots()
		ax.set_title('$c^2$')
		ax.set_xlabel(xname)
		for i,fld in enumerate(fld_list):
			ax.plot(r[i],fld.c2,linestyle)
			
	if q=='hor':	
		fig,ax = subplots()
		ax.set_title('$h/r$')
		ax.set_xlabel(xname)
		for i,fld in enumerate(fld_list):
			ax.plot(r[i],fld.hor,linestyle)
		
	if q=='E':
		fig,(ax1,ax2,ax3,ax4)=subplots(4,sharex=True)
		ax1.set_title('Eccentricity')
		ax1.set_ylabel('$e_x$')
		ax2.set_ylabel('$e_y$')
		ax3.set_ylabel('e')
		ax4.set_ylabel('$\\omega/\\pi$')
		ax4.set_xlabel(xname)
		
		for i,fld in enumerate(fld_list):
			ax1.plot(r[i],real(fld.E),linestyle)
			ax2.plot(r[i],imag(fld.E),linestyle)
			ax3.plot(r[i],abs(fld.E),linestyle)
			ax4.plot(r[i],angle(fld.E)/pi,linestyle)
	
	if q=='e':
		fig,ax = subplots()
		ax.set_title('$e$')
		ax.set_xlabel(xname)
		for i,fld in enumerate(fld_list):
			ax.plot(r[i],abs(fld.E),linestyle)
			
	if q=='w':
		fig,ax = subplots()
		ax.set_title('$\\omega/\\pi$')
		ax.set_xlabel(xname)
		for i,fld in enumerate(fld_list):
			ax.plot(r[i],angle(fld.E)/pi,linestyle)
			
			
	if q=='ex':
		fig,ax = subplots()
		ax.set_title('$e_x$')
		ax.set_xlabel(xname)	
		for i,fld in enumerate(fld_list):
			ax.plot(r[i],real(fld.E),linestyle)
		
	if q=='ey':
		fig,ax = subplots()
		ax.set_title('$e_y$')
		ax.set_xlabel(xname)	
		for i,fld in enumerate(fld_list):
			ax.plot(r[i],imag(fld.E),linestyle)
			

	legend(leg_str,loc='best')
	show()
	

class star():
	def __init__(self,r,dbar,m,soft):
			self.rs = 2*pi * trapz(dbar*r*r,x=r)
			phi = linspace(0,pi,2000)
			self.gr = zeros(r.shape)
			self.gp = zeros(r.shape,dtype='complex')
			self.phi = zeros(r.shape)
			for i in range(len(r)):
				rsoft = r[i]*r[i] + soft*soft
				q = self.rs / sqrt(rsoft)	
				self.phi[i] = -2/sqrt(rsoft) * trapz(cos(m*phi)/sqrt(1+q*q+2*q*cos(phi)),x=phi)
				self.gr[i] = -2*self.rs*r[i] / rsoft**2 * trapz(cos(m*phi)*(cos(phi) +q)*(1+q*q+2*q*cos(phi))**(-1.5),x=phi)
				self.gp[i] = -1j*m * self.phi[i] / sqrt(rsoft)

def E_pred(r,Ei,Eo,beta,alpha_b,alpha_s,gamma_s,gamma_b,bc,ogilvie=False):
# bc = 0 -> E(ri) = Ei & E(ro) = Eo
# bc = 1 -> E'(ri) = 0 & E(ro) = Eo
# bc = 2 -> E(ri) = Ei & E'(ro) = 0
# bc = 3 -> E(ro) = Eo & E'(ro) = 0 
# bc = 4 -> E(ri) = Ei & E'(ri) = 0 


	ri = r[0]
	ro = r[-1]
	chi = ri/ro
#	eta = (1 + 1j*nu)/(1 + nu*nu)
#	gam = sqrt( 1 + beta*beta/4 + beta*(1 - eta))


#	ai = -1 - beta/2 - gam
#	ao = -1 - beta/2 + gam

	gamma_s += beta
	gamma_b += beta
	if ogilvie:
		gam = sqrt( 1 + beta*beta/4)
		ai = -1 - beta/2 - gam
		ao = -1 - beta/2 + gam		
	else:
	
		a2= 1-1j*alpha_b
	
#		a0 = 2*beta-1j*alpha_s*(2-2.5*gamma_s+1.5*beta)
#		a1 = beta+3-1j*alpha_b*(gamma_b+1.5)-1j*alpha_s*(3*(gamma_s-beta) - .5)
	
		a0 = 2*beta-1j*.5*alpha_s*(3 + 3*beta - gamma_s)
		a1 = beta+3-1j*alpha_b*(gamma_b+1.5)-1j*3*alpha_s*(gamma_s-beta)
		
		a0 /= a2
		a1 /= a2
	
		ai = .5*(1-a1 - sqrt((a1-1)**2 -4*a0))
		ao = .5*(1-a1 + sqrt((a1-1)**2 -4*a0))

	print 'Density power law is: ',beta
	print 'Inner power law is: ', ai
	print 'Outer power law is: ', ao
	
# 	if bc not in ['Dirichlet','Neumann','dirichlet','neumann']:
# 		print 'Not a valid boundary condition'
# 		return
		
	if bc not in [0,1,2,3,4]:
		print 'Not a valid b.c, using bc=0 as default'
		bc = 0
		
		
# 	if bc == 0:
# 		A = (Ei - Eo*pow(chi,-ao))/(1 - pow(chi,ai-ao))
# 		B = (Eo - Ei*pow(chi,ai))/(1- pow(chi,ai-ao))
# 		
# 	if bc==1:
# 		A = -Eo*pow(chi,1-ao)/(1 - pow(chi,1+ai-ao))
# 		B = Eo/(1 - pow(chi,1+ai-ao))
# 	
# 	if bc==2:
# 		A = Ei/(1-pow(chi,ai-ao-1))
# 		B = -Ei*pow(chi,ai-1)/(1- pow(chi,ai-ao-1))
		
	
	if bc == 0:
		A = (Ei - Eo*pow(chi,ao))/(1 - pow(chi,ao-ai))
		B = (Eo - Ei*pow(chi,-ai))/(1- pow(chi,ao-ai))
		
	if bc==1:
		A = Eo*ao*pow(chi,ao)/(ao*pow(chi,ao-ai)-ai)
		B = -Eo*ai/(ao*pow(chi,ao-ai)-ai)
	
	if bc==2:
		A = Ei*ao/(ao - ai*pow(chi,ao-ai))
		B = -Ei*ai*pow(chi,-ai)/(ao - ai*pow(chi,ao-ai))
	
	if bc==3:
		A = Eo*ao*pow(chi,-ai)/(ao-ai)
		B = -Eo*ai/(ao-ai)
	if bc==4:
		A = Ei*ao/(ao-ai)
		B = -Ei*ai*pow(chi,-ao)/(ao-ai)
		
		
	print 'Inner coefficient: ', A
	print 'Outer coefficient: ', B
	
		
	out = A*pow(r/ri,ai) + B*pow(r/ro,ao)	
#	out *= ((((abs(real(out))<1e-10) & (abs(imag(out))<1e-10)).astype(int) + 1 ) % 2 )
	return out
	
	
# 	if bc in ['Dirichlet','dirichlet']:
# 		faci = Ei * pow(r/ri,plaw) * ( ( 1 - pow(r/ro,-gam))/(1 - pow(ri/ro,-gam)))
# 		faco = Eo * pow(r/ro,plaw) * ( ( 1 - pow(r/ri,-gam))/(1 - pow(ro/ri,-gam)))
# 		return faci + faco
# 	
# 	if bc in ['Neumann','neumann']:
# 		bfrac = beta/gam
# 		norm = Eo * pow(r/ro,plaw)
# 		fac1 = ( pow(r/ri,-gam) * ( 1 - bfrac) + 1 + bfrac)
# 		fac2 = (pow(ro/ri,-gam) * ( 1 - bfrac) + 1 + bfrac)
# 		return norm * fac1 / fac2
	
# 	faco = Eo*pow(ro,.5*(beta + bfac))*pow(ri,bfac)*(1-pow(r/ri,bfac))
# 	faci =  Ei*pow(ri,.5*(beta + bfac))*pow(ro,bfac)*(1-pow(r/ro,bfac))
# 	norm = pow(r,-1-.5*(beta + bfac))
# 	norm /= (pow(ri,bfac) - pow(ro,bfac))
	
	
	
	
def eccen_plots(r,bvals,Ei,Eo,alpha_b,alpha_s,gamma_s,gamma_b,bc,linestyle='-',logscale=True,yscale=False,ogilvie=False):
	
	epred = [E_pred(r,Ei,Eo,b,alpha_b,alpha_s,gamma_s,gamma_b,bc,ogilvie=ogilvie) for b in bvals]

		
	nu_str = '$\\alpha$ = %.1e' % alpha_s
	nu_str += ' and $d\ln\\nu/d\lnr$ = %.2f' % gamma_s
	if bc not in [0,1,2,3,4]:
		bc = 0
	
	if bc==0:
		tstr = '$E(r_i)=$'+str(Ei)+' and $E(r_o)=$'+str(Eo) + ', '+ nu_str
	if bc==1:
		tstr = '$d_r E(r_i)= 0 $ and $E(r_o)=$'+str(Eo)+ ', '+ nu_str
	if bc==2:
		tstr = '$E(r_i)=$'+str(Ei)+' and $d_rE(r_o)=0$'+ ', '+ nu_str
	if bc==3:
		tstr = '$E(r_o)=$'+str(Eo)+' and $d_rE(r_o)=0$'+ ', '+ nu_str
	if bc==4:
		tstr = '$E(r_i)=$'+str(Ei)+' and $d_rE(r_i)=0$'+ ', '+ nu_str
	
	fig,(ax1,ax2)=subplots(2,sharex=True)
	ax1.set_title(tstr)
	ax2.set_ylabel('$\\omega/ \\pi$')
	
	if logscale:
		if bc == 2 or bc==4:
			ax1.set_ylabel('$\ln (e/e_i)$')
		else:
			ax1.set_ylabel('$\ln (e/e_o)$')
		ax2.set_xlabel('$\ln r$')
	else:
		if bc==2 or bc==4:
			ax1.set_ylabel('$e/e_i$')
		else:
			ax1.set_ylabel('$e/e_o$')
		ax2.set_xlabel('$r$')
	
	for i,ep in enumerate(epred):
		if bc==2 or bc==4:
			dat = abs(ep)/abs(ep[0])
		else:
			dat = abs(ep)/abs(ep[-1])
			
		if logscale:
				
			ax1.plot(log(r),log(dat),linestyle,label='$\\beta=$'+str(bvals[i]))
		else:
			ax1.plot(r,dat,linestyle,label='$\\beta=$'+str(bvals[i]))
	
	
	ax1.legend(loc='best')
#	gca().set_color_cycle(None)	
	for i,ep in enumerate(epred):
		if logscale:
			ax2.plot(log(r),angle(ep)/pi,linestyle)
		else:
			ax2.plot(r,angle(ep)/pi,linestyle)
		
#	ax2.set_ylim(-1.1,1.1)
	
	if yscale:
		if logscale:
			ax1.set_ylim(-6,2)
	show()
	
	
def find_surface_density( disk_mass, rin, rout, beta):
	
	if beta == -2:
		fac = 2 * pi * log(rout/rin)
	else:
		fac = 2*pi * (rout**(beta+2) - rin**(beta+2))
		fac /= (beta + 2)
	
	return disk_mass / fac
	
	
def poisson_kernel(r,hor,eps):
	
	H = hor*r;
	phi = linspace(-pi,pi,200)
	Nr = len(r)
	
	kernel = zeros((Nr,Nr))
	
	for i in range(Nr):
		for j in range(Nr):
			chi  = (r[i]*r[i] + r[j]*r[j] - 2*r[i]*r[j]*cos(phi) + eps*eps*H[j]*H[j])/(4*H[j]*H[j])
			kernel[i,j] = 	trapz((1.0/sqrt(2*pi*H[j]*H[j]))*kn(0,chi)*exp(chi),x=phi)
	
	figure()
	imshow(kernel,center='origin')
	colorbar()
	
	return kernel

def compare_eccen_pred(r,edat, epred,epredog,tstr,rsample=200,logr=False,norm=False,secondtstr=None):
	
	if logr:
		xname = '$\ln r$'
	else:
		xname = '$r$'
	
	if norm:
		epnorm = abs(epred).max()
		eponorm = abs(epredog).max()
		ednorm = abs(edat).max()
		
	else:
		epnorm = 1
		eponorm = 1
		ednorm = 1

	fig,((ax1,ax3),(ax2,ax4)) = subplots(2,2,sharex='col')
		
#	fig,(ax1,ax2)=subplots(2,sharex=True)
	ax1.set_title(tstr,fontsize='large')
	
	if secondtstr != None:
		ax3.set_title(secondtstr,fontsize='large')
		
	if norm:
		ax1.set_ylabel('$e/e_{max}$',fontsize='large')
	else:
		ax1.set_ylabel('$e$',fontsize='large')
		
	ax2.set_ylabel('$\\omega / \\pi$',fontsize='large')
	ax2.set_xlabel(xname,fontsize='large')
	
	
	ax4.set_xlabel(xname,fontsize='large')
	
	
	ax3.set_ylabel('$e_y$',fontsize='large')
	ax4.set_ylabel('$e_x$',fontsize='large')	

	ax1.plot(r,abs(epred)/epnorm,'b-',label=r'Prediction')
	ax1.plot(r,abs(epredog)/eponorm,'r-',label=r'Ogilvie Prediction')
	ax1.plot(r[::rsample],abs(edat[::rsample])/ednorm,'kx',label=r'Data')
	ax2.plot(r,angle(epred)/pi,'b-',label=r'Prediction')
	ax2.plot(r,angle(epredog)/pi,'r-',label=r'Ogilvie Prediction')
	ax2.plot(r[::rsample],angle(edat[::rsample])/pi,'kx',label=r'Data')
	
	
	
	
	ax3.plot(r,imag(epred),'b-',label=r'Prediction')
	ax3.plot(r,imag(epredog),'r-',label=r'Ogilvie Prediction')
	ax3.plot(r[::rsample],imag(edat[::rsample]),'kx',label=r'Data')
	
	ax4.plot(r,real(epred),'b-',label=r'Prediction')
	ax4.plot(r,real(epredog),'r-',label=r'Ogilvie Prediction')
	ax4.plot(r[::rsample],real(edat[::rsample]),'kx',label=r'Data')
	
	ax4.legend(loc='best')
	
	return 

def prediction_results(bc,ei=0,eo=.1,normalize=False):
	beta_vals = [-1.5, -.75, -.5, 0]
	salpha_vals = [ .03, .003, .0003,.00003]
	balpha_vals = [ -.02, -.002, -.0002,-.00002]
	
	
# 	if bc in [2,4]:
# 		ei = eo
# 		eo = 0
# 		
	
		
	fld=[]
	beta=[]
	alpha=[]
	epred=[]
	epredog=[]
	eresids=[]
	wresids=[]
	
	
	
	if bc==0:
		bctstr = '$E(r_i)=$'+str(ei)+' and $E(r_o)=$'+str(eo) 
	if bc==1:
		bctstr = '$d_r E(r_i)= 0 $ and $E(r_o)=$'+str(eo)
	if bc==2:
		bctstr = '$E(r_i)=$'+str(ei)+' and $d_rE(r_o)=0$'
	if bc==3:
		bctstr = '$E(r_o)=$'+str(eo)+' and $d_rE(r_o)=0$'
	if bc==4:
		bctstr = '$E(r_i)=$'+str(ei)+' and $d_rE(r_i)=0$'
	
	
	
	for i in range(len(beta_vals)):
		for j in range(len(salpha_vals)):
			b = beta_vals[i]
			a = salpha_vals[j]
			beta.append(b)
			alpha.append(a)
			
			
			dir_name = 'analytic_'+str(b)+'_'+str(a)+'/'
			
			print dir_name
			tfld=Field(1,outdir=dir_name)
			fld.append(tfld)
			tepred = E_pred(tfld.nlr,ei,eo,b,balpha_vals[j],a,tfld.ind_nu,tfld.ind_nu,bc)
			tepredog = E_pred(tfld.nlr,ei,eo,b,balpha_vals[j],a,tfld.ind_nu,tfld.ind_nu,bc,ogilvie=True)
			epred.append(tepred)
			epredog.append(tepredog)
			
			tstr = '$\\beta = %.2f$' % b
			tstr += ', $\\alpha = %.2e$' % a
						
			eresids.append(l2norm(abs(tepred)/abs(tepred).max() - abs(tfld.E)/abs(tfld.E).max()))
			wresids.append(l2norm(angle(tepred)/pi - angle(tfld.E))/pi)
			
			compare_eccen_pred(tfld.r,tfld.E,tepred,tepredog,tstr,rsample=100,norm=normalize,logr=True,secondtstr=bctstr)
	
		
	return fld, epred, beta, alpha,eresids,wresids		
			

def fit_growth_rates(tstart=None,tend=None,outdir=''):
	dat = loadtxt(outdir + 'history.dat')
	fig, ((ax1,ax2),(ax3,ax4)) = subplots(2,2)
	
	
	if tend != None:
		dat=dat[dat[:,0]<=tend,:]
	if tstart != None:
		dat = dat[dat[:,0]>=tstart,:]
		
	
	evals = dat[:,7]
	
	t = dat[:,0]
	
	norm = len(t)
	
	freq = 2*pi*fft.rfftfreq(len(dat[:,0]),diff(dat[:,0]).mean())
	
	s = polyfit(t,log(evals),1)
	
	growth_rate = s[0]
	
	expfit = exp( s[1] + t*s[0])
	
	
	ft = fft.rfft(evals/expfit) 
	
	ax3.set_title('Residual Power Spectrum',fontsize='large')
	ax3.set_xlabel('$\log_{10} \\Omega_p$',fontsize='large')
	ax3.set_ylabel('FT power',fontsize='large')
	ax3.plot(log10(freq),abs((ft/norm)**2),'-x')
	
	
	
	ind = abs(ft**2) == abs(ft**2)[1:].max()
	
	omega = freq[ind][0]
	
	
	for i,x in enumerate(ft):
		if i != 0 and freq[i] != omega:
			ft[i] = 0	
			
	efit = fft.irfft(ft,len(t))


	
	
	efit_final= efit*expfit
	
	
	
	
	
	
	
	ax1.set_title('Total Fit',fontsize='large')
	ax1.set_xlabel('t',fontsize='large')
	ax1.set_ylabel('e',fontsize='large')
	ax1.plot(t,evals)
	ax1.plot(t,expfit,'--',label=r'$\gamma$ = %.2e' % growth_rate)
	ax1.plot(t,efit_final,label=r'$\Omega_p$ = %.2e' % omega )
	ax1.legend(loc='best')
	
	
	ax2.set_title('One Frequency Pattern Speed Fit',fontsize='large')
	ax2.set_xlabel('t',fontsize='large')
	ax2.set_ylabel('e/ exponential fit',fontsize='large')
	ax2.plot(t,evals/expfit)
	ax2.plot(t,efit,label=r'$\Omega_p$ = %.2e' % omega)
	
	
	
	
	
	ax4.set_title('Exponential fit',fontsize='large')
	ax4.set_xlabel('t',fontsize='large')
	ax4.set_ylabel('e',fontsize='large')
	ax4.plot(t,evals)
	ax4.plot(t,expfit,label=r'$d\ln e/ dt$ = %.2e' % growth_rate)
	
	
	
	
	
	fig2, ((ax12,ax22),(ax32,ax42)) = subplots(2,2)
	
	
	
	if tstart == None:
		rs = dat[1:,4]
		t=t[1:]
	else:
		rs = dat[:,4]
		t=t[:]
		
		
	
	s_star = polyfit(t,log(rs),1)
	
	growth_rate_star = s_star[0]
	
	expfit_star = exp( s_star[1] + t*s_star[0])
	
	
	ft_star = fft.rfft(rs/expfit_star) 
	
	ax32.set_title('Residual Power Spectrum',fontsize='large')
	ax32.set_xlabel('$\log_{10} \\Omega_\star$',fontsize='large')
	ax32.set_ylabel('FT power',fontsize='large')
	ax32.plot(log10(freq),abs((ft_star/norm)**2),'-x')


	maxft = sort(abs(ft_star**2)[1:])[::-1][:2]

		
	ind_star_1 = abs(ft_star**2) == maxft[0]
	ind_star_2 = abs(ft_star**2) == maxft[1]

	omega_star_1 = freq[ind_star_1][0]
	omega_star_2 = freq[ind_star_2][0]
	
	omstr_star = 'low $\Omega_\star$ = %.2e' % min([omega_star_1,omega_star_2])
	omstr_star += ', high $\Omega_\star$ = %.2e' % max([omega_star_1,omega_star_2])
	
	print omega_star_1,omega_star_2
	
	for i,x in enumerate(ft_star):
		if i != 0 and freq[i] != omega_star_1 and freq[i] != omega_star_2:
			ft_star[i] = 0	
	

	efit_star = fft.irfft(ft_star,len(t))
	
	
	efit_final_star = efit_star*expfit_star
	
	
	
	
	
	
	
	ax12.set_title('Total Fit',fontsize='large')
	ax12.set_xlabel('t',fontsize='large')
	ax12.set_ylabel('$r_\star$',fontsize='large')
	ax12.plot(t,rs)
	ax12.plot(t,expfit_star,'--',label=r'$\gamma_\star$ = %.2e' % growth_rate_star)
	ax12.plot(t,efit_final_star,label=omstr_star )
	ax12.legend(loc='best')
	
	
	ax22.set_title('Two Frequency Pattern Speed Fit',fontsize='large')
	ax22.set_xlabel('t',fontsize='large')
	ax22.set_ylabel('$r_\star$ / exponential fit',fontsize='large')
	ax22.plot(t,rs/expfit_star)
	ax22.plot(t,efit_star,label=omstr_star)
	
	
	
	
	
	ax42.set_title('Exponential fit',fontsize='large')
	ax42.set_xlabel('t',fontsize='large')
	ax42.set_ylabel('$r_\star$',fontsize='large')
	ax42.plot(t,rs)
	ax42.plot(t,expfit_star,label=r'$d\ln e/ dt$ = %.2e' % growth_rate_star)
	
	
	
	return (growth_rate, omega),(growth_rate_star,omega_star_1,omega_star_2)
	
	
	
	
	

		