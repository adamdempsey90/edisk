class Field():
	def __init__(self,t,outdir='',loadrhs=False,rhsnum=0,NG=1):	
		disk = loadtxt(outdir+'disk.dat')
		self.r = disk[:,0]
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
		self.dbar = dat[:,9]
		self.E = (2*self.v - 1j*self.u) / (2*self.omk*self.r)
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
	def plot(self,q,linestyle='-'):
		if q=='u':
			fig,(ax1,ax2)=subplots(2,sharex=True)
			ax1.set_title('u')
			ax1.plot(self.r,real(self.u),linestyle,label=r'$Re(u)$')
			ax1.set_ylabel('Re(u)')
			ax2.plot(self.r,imag(self.u),linestyle,label=r'$Im(u)$')
			ax2.set_ylabel('Im(u)')
			ax2.set_xlabel('r')
		if q=='v':
			fig,(ax1,ax2)=subplots(2,sharex=True)
			ax1.set_title('v')
			ax1.plot(self.r,real(self.v),linestyle,label='$Re(v)$')
			ax1.set_ylabel('Re(v)')
			ax2.plot(self.r,imag(self.v),linestyle,label='$Im(v))$')
			ax2.set_ylabel('Im(v)')
			ax2.set_xlabel('r')
		if q=='sig':
			fig,(ax1,ax2)=subplots(2,sharex=True)
			ax1.set_title('$\\sigma / <\\Sigma>$')
			ax1.plot(self.r,real(self.sig),linestyle,label='$Re(\\sigma))$')
			ax1.set_ylabel('Re($\\sigma$)')
			ax2.plot(self.r,imag(self.sig),linestyle,label='$Im(\\sigma))$')
			ax2.set_ylabel('Im($\\sigma$)')
			ax2.set_xlabel('r')
		if q=='vybar':
			fig,ax = subplots()
			ax.set_title('$<v_\\phi>$')
			ax.plot(self.r,self.vyb,linestyle)
			ax.set_xlabel('r')
			
		if q=='dbar':
			fig,ax = subplots()
			ax.set_title('$<\\Sigma>$')
			ax.plot(self.r,self.dbar,linestyle)
			ax.set_xlabel('r')
			
		if q=='omk':
			fig,ax = subplots()
			ax.set_title('$\\Omega_k$')
			ax.plot(self.r,self.omk,linestyle)
			ax.set_xlabel('r')
		if q=='nu':
			fig,ax = subplots()
			ax.set_title('$\\nu$')
			ax.plot(self.r,self.nu,linestyle)
			ax.set_xlabel('r')
		if q=='c2':
			fig,ax = subplots()
			ax.set_title('$c^2$')
			ax.plot(self.r,self.c2,linestyle)
			ax.set_xlabel('r')
		if q=='hor':	
			fig,ax = subplots()
			ax.set_title('$h/r$')
			ax.plot(self.r,self.hor,linestyle)
			ax.set_xlabel('r')
			
		if q=='E':
			fig,(ax1,ax2,ax3,ax4)=subplots(4,sharex=True)
			ax1.plot(self.r,real(self.E),linestyle)
			ax1.set_ylabel('$e_x$')
			ax2.plot(self.r,imag(self.E),linestyle)
			ax2.set_ylabel('$e_y$')
			ax3.plot(self.r,abs(self.E),linestyle)
			ax3.set_ylabel('e')
			ax4.plot(self.r,angle(self.E),linestyle)
			ax4.set_ylabel('$\\omega$')
			ax4.set_xlabel('r')
		
		if q=='dtu':
			fig,(ax1,ax2)=subplots(2,sharex=True)
			ax1.set_title('dtu')
			ax1.plot(self.r[self.NG:-self.NG],real(self.dtu),linestyle,label=r'$Re(dtu)$')
			ax1.set_ylabel('Re(dtu)')
			ax2.plot(self.r[self.NG:-self.NG],imag(self.dtu),linestyle,label=r'$Im(dtu)$')
			ax2.set_ylabel('Im(dtu)')
			ax2.set_xlabel('r')
		if q=='dtv':
			fig,(ax1,ax2)=subplots(2,sharex=True)
			ax1.set_title('dtv')
			ax1.plot(self.r[self.NG:-self.NG],real(self.dtv),linestyle,label='$Re(dtv)$')
			ax1.set_ylabel('Re(dtv)')
			ax2.plot(self.r[self.NG:-self.NG],imag(self.dtv),linestyle,label='$Im(dtv))$')
			ax2.set_ylabel('Im(dtv)')
			ax2.set_xlabel('r')
		if q=='dts':
			fig,(ax1,ax2)=subplots(2,sharex=True)
			ax1.set_title('dts')
			ax1.plot(self.r[self.NG:-self.NG],real(self.dts),linestyle,label='$Re(dts)$')
			ax1.set_ylabel('Re(dts)')
			ax2.plot(self.r[self.NG:-self.NG],imag(self.dts),linestyle,label='$Im(dts)$')
			ax2.set_ylabel('Im(dts)')	
		
		show()
		
		
		
def animate(q,t,dt=1,linestyle='-'):
	fld0 = Field(0)
	dat = zeros((fld0.nr,len(t)))
	for i,j in enumerate(t):
		print 'Loading time ',j
		fld = Field(j)
		if q=='u':
			dat[:,i] = fld.u
		if q=='v':
			dat[:,i] = fld.u
		if q=='sig':
			dat[:,i] = fld.u
		if q=='E':
			dat[:,i] = fld.E
		
	
	if q=='u':
		fig,(ax1,ax2)=subplots(2,sharex=True)
		ax1.set_title('u')
		ax1.set_ylabel('Re(u)')
		ax2.set_ylabel('Im(u)')
		ax2.set_xlabel('r')
		
		for i in range(len(t)):
			fig.clf()
			ax1.plot(fld0.r,real(dat[:,i]),linestyle)
			ax2.plot(fld0.r,imag(dat[:,i]),linestyle)
			fig.canvas.draw()
			
		
	if q=='v':
		fig,(ax1,ax2)=subplots(2,sharex=True)
		ax1.set_title('v')
		ax1.set_ylabel('Re(v)')
		ax2.set_ylabel('Im(v)')
		ax2.set_xlabel('r')
		for i in range(len(t)):
			fig.clf()
			ax1.plot(fld0.r,real(dat[:,i]),linestyle)
			ax2.plot(fld0.r,imag(dat[:,i]),linestyle)
			fig.canvas.draw()
	if q=='sig':
		fig,(ax1,ax2)=subplots(2,sharex=True)
		ax1.set_title('$\\sigma / <\\Sigma>$')
		ax1.set_ylabel('Re($\\sigma$)')
		ax2.set_ylabel('Im($\\sigma$)')
		ax2.set_xlabel('r')
		for i in range(len(t)):
			fig.clf()
			ax1.plot(fld0.r,real(dat[:,i]),linestyle)
			ax2.plot(fld0.r,imag(dat[:,i]),linestyle)
			fig.canvas.draw()
	if q=='E':
		fig,(ax1,ax2,ax3,ax4)=subplots(4,sharex=True)
		ax1.set_ylabel('$e_x$')
		ax2.set_ylabel('$e_y$')
		ax3.set_ylabel('e')
		ax4.set_ylabel('$\\omega$')
		ax4.set_xlabel('r')	
		for i in range(len(t)):
			fig.clf()
			ax1.plot(fld0.r,real(dat[:,i]),linestyle)
			ax2.plot(fld0.r,imag(dat[:,i]),linestyle)
			ax3.plot(fld0.r,abs(dat[:,i]),linestyle)
			ax4.plot(fld0.r,angle(dat[:,i]),linestyle)
			fig.canvas.draw()
	
			