class Field():
	def __init__(self,t,outdir='',loadrhs=False,rhsnum=0,NG=1,etol=1e-8):	
		disk = loadtxt(outdir+'disk.dat')
		self.r = disk[:,0]
		self.nlr = 10**self.r
		self.nr = len(self.r)
		self.dr = diff(self.r)[0]
		self.width = self.nr*self.dr
		self.hor = disk[:,1]
		self.c2 = disk[:,2]
		self.nu = disk[:,3]
		dat = loadtxt(outdir + 'output_'+str(t)+'.dat')
		self.u = dat[:,1]+1j*dat[:,2]
		self.v = dat[:,3] + 1j*dat[:,4]
		self.sig = dat[:,5] + 1j*dat[:,6]
		self.vyb = dat[:,7]
		self.omk = dat[:,8]
		self.omk0 = pow(self.nlr,-1.5)
		self.dbar = dat[:,9]
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
			self.dtu = rhs[:,1]+rhs[:,2]*1j
			self.dtv = rhs[:,3] + rhs[:,4]*1j
			self.dts = rhs[:,5] + rhs[:,6]*1j
		else:
			self.dtu = 0
			self.dtv = 0
			self.dts = 0
	def plot(self,q,linestyle='-',logr=True):
		if logr:
			r = self.r
			xname = '$\ln r$'
		else:
			r = 10**self.r
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
			ax4.plot(r,angle(self.E),linestyle)
			ax4.set_ylabel('$\\omega$')
			ax4.set_xlabel(xname)
		
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
		
def animate(q,t,dt=1,linestyle='-',dat=None,fld0=None,logr=True):
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
		
	
	if logr:
		r = fld0.r
		xstr = '$\ln r$'
	else:
		r = fld0.nlr
		xstr = 'r'
		
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
		ax4.set_ylabel('$\\omega$')
		ax4.set_xlabel(xstr)	
		
		l1,= ax1.plot(r,real(fld0.E),linestyle)
		l2,= ax2.plot(r,imag(fld0.E),linestyle)
		l3,= ax3.plot(r,abs(fld0.E),linestyle)
		l4,= ax4.plot(r,angle(fld0.E),linestyle)
		
		l1range = (real(dat).min(),real(dat).max())
		l2range = (imag(dat).min(),imag(dat).max())
		l3range = (abs(dat).min(),abs(dat).max())
		l4range = (angle(dat).min(),angle(dat).max())
		fig.canvas.draw()
		for i in range(len(t)):
			ax1.set_title('Eccentricity, t='+str(t[i]*dt))
			l1.set_ydata(real(dat[:,i]))
			l2.set_ydata(imag(dat[:,i]))
			l3.set_ydata(abs(dat[:,i]))
			l4.set_ydata(angle(dat[:,i]))

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
	