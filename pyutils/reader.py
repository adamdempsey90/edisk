class Field():
	def __init__(self,t,outdir='',loadrhs=False,rhsnum=0,NG=1,etol=1e-8):	
		disk = loadtxt(outdir+'disk.dat')
		self.r = disk[:,0]
		self.nlr = disk[:,1]
		self.nr = len(self.r)
		self.dr = diff(self.r)[0]
		self.width = self.nr*self.dr
		self.hor = disk[:,2]
		self.c2 = disk[:,3]
		self.nu = disk[:,4]
		dat = loadtxt(outdir + 'output_'+str(t)+'.dat')
		self.u = dat[:,2]+1j*dat[:,3]
		self.v = dat[:,4] + 1j*dat[:,5]
		self.sig = dat[:,6] + 1j*dat[:,7]
		self.vyb = dat[:,8]
		self.omk = dat[:,9]
		self.omk0 = pow(self.nlr,-1.5)
		self.dbar = dat[:,10]
		self.E = (2*self.v - 1j*self.u) / (2*self.vyb)
		

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
		else:
			self.dtu = 0
			self.dtv = 0
			self.dts = 0
	def plot(self,q,linestyle='-',logr=True):
	
		if q not in ['u','v','sig','E','nu','c2','hor','omk','dbar','vybar','dtu','dtv','dts','e','w','ex','ey']:
			print 'Not Valid Variable Name'
			return
			
		if logr:
			r = self.r
			xname = '$\ln r$'
		else:
			r = self.nlr
			xname = '$r$'
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
			fig,ax = subplots()
			ax.set_title('$\\nu$')
			ax.plot(r,self.nu,linestyle)
			ax.set_xlabel(xname)
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
			ax.plot(r,abs(self.E),linestyle)
			ax.set_xlabel(xname)
		if q=='w':
			fig,ax = subplots()
			ax.set_title('$e$')
			ax.plot(r,angle(self.E),linestyle)
			ax.set_xlabel(xname)
		if q=='ex':
			fig,ax = subplots()
			ax.set_title('$e$')
			ax.plot(r,real(self.E),linestyle)
			ax.set_xlabel(xname)	
		if q=='ey':
			fig,ax = subplots()
			ax.set_title('$e$')
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
		
		
		show()
		
		
		
		
	def draw_ellipse(self,num_ellipse,(xc,yc)=(0,0),Nph=200):
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
		for i in range(self.nr)[::self.nr/num_ellipse]:
			plot(x[:,i],y[:,i],'-k')
			
		return x,y
	
	def plotdens(self,Nph=200,starpos=(0,0)):
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
				dat[:,i] = abs(fld.E)
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
	

def animate_ellipse(t,num_ellipse,(xc,yc)=(0,0),Nph=200):
		
		
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


def animate_real(q,t,xlims=None,ylims=None,Nph=200):

	if q not in ['dens','vx','vy','E']:
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
	
	for i,j in enumerate(t):
		print 'Loading t = ' + str(i)
		fld=Field(j)
		for p in range(Nph):
			if q=='vx':
				dat[:,p,i] =  2*real(fld.u*exp(-1j*phi[p]))
			if q=='vyp':
				dat[:,p,i] =   2*real(fld.v*exp(-1j*phi[p]))
			if q=='vy':
				dat[:,p,i] =   2*real(fld.v*exp(-1j*phi[p])) + fld.vyb
			if q=='dens':
				dat[:,p,i] =  fld.dbar*(1 + 2*real(fld.sig*exp(-1j*phi[p])))/fld.dbar - 1
			if q=='E':
				dat[:,p,i] =  2*abs(real(fld.E*exp(-1j*phi[p])))
	
	fig=figure()
	for i,j in enumerate(t[1:]):
		fig.clear()
		pcolormesh(x,y,dat[:,:,i])
		colorbar()
		title(q + ',    t = ' + str(j))
		fig.canvas.draw()
		
	
	
	return dat,x,y,fld0


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

