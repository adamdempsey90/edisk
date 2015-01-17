class fld():
	def __init__(self,m,lr,Mdisk,h,beta,alphas,alphab,flare):
		self.m = m
		self.lr = lr
		self.Nr = len(lr)
		self.dr = diff(lr)[0]
		self.r = exp(lr)
		self.E = zeros((self.Nr,1),dtype='complex')
		self.sig = zeros((self.Nr,1),dtype='complex')
		
		self.omk = pow(self.r,-1.5)
		self.dbar = pow(self.r,beta)
		self.dbar *= Mdisk/(2*pi*trapz(self.dbar*self.r,x=self.r))
		self.hor = h*pow(self.r,flare)
		self.nus = alphas * self.hor * self.hor * self.omk * self.r * self.r
		self.nub = alphab * self.hor * self.hor * self.omk * self.r * self.r
		self.gamma = 2*(1+flare) + q + beta
		self.alpha_s = alphas
		self.alpha_b = alphab
	
	def steady_state(self,bc):
		G = [zeros((2,1),dtype='complex')] * self.Nr
		H = [zeros((2,2),dtype='complex')] * self.Nr
		UK = zeros((2,3),dtype='complex')
		HG = zeros((2,3),dtype='complex')

# Forward Solve
		for i in range(self.Nr):
			
			M,U,L,K = load_matrices(i,self,bc):
		
			if i!=0:
				K -= dot(L,G[i-1])		
		
			UK[:,:2] = -U
			UK[:,-1] = K
			
			if i!=0:
				M += dot(L,H[i-1])
			
		 	temp = solve(M,UK)
			H[i] = temp[:,:2]
			G[i] = temp[:,-1]

# Backward Substitution		
		
		
		
		self.E[-1] = G[-1][0]
		self.sig[-1] = G[-1][0]
		
		for i in range(self.Nr-1)[::-1]:
			G[i] += dot(H,G[i+1])
			
			self.E[i] = G[i][0]
			self.sig[i] = G[i][0]
			
		return		
		
	def load_matrices(i,self,bc):
		
		A = zeros((2,2),dtype='complex')
		B = zeros((2,2),dtype='complex')
		C = zeros((2,2),dtype='complex')
		M = zeros((2,2),dtype='complex')
		U = zeros((2,2),dtype='complex')
		L = zeros((2,2),dtype='complex')
		K = zeros((2,1),dtype='complex')
		
		r = self.r[i]
		omk =self.omk[i]
		nus = self.nus[i]
		nub = self.nub[i]
		h2o2 = self.hor[i]**2 r**2 * omk**2 
		
		
		
		if i==0:
			M[0][0] = 1
			M[1][1] = 1
			if bc[0]=='fixed':
				K[0] =  bc[1]
				
			if bc[0]=='innner':
				U[0][0] = -1
				U[1][1] = -1
			if bc[0]=='outer':
				K[0] = bc[1]
			
		
		
		
		elif i==self.Nr-1:
			M[0][0] = 1
			M[1][1] = 1
			if bc[0]=='fixed':
				K[0] =  bc[2]
				
			if bc[0]=='innner':
				L[0][0] = -1
				L[1][1] = -1
			if bc[0]=='outer':
				K[0] = bc[2]
		
		
		else:
			A[0][0] = -nus*(self.gamma - 3) * .5 * (omk/r)
			A[0][1] = (4*1j++3*self.alpha_s)*h2o2/(2*r) 
			A[1][0] = -1j*self.beta*omk
			A[1][1] = 1j*omk
			
			B[0][0] = .5*(3*(nus + nub) + 2*(nub + 3*nubs)*self.gamma)*omk
			B[0][1] = (1j-3*self.alpha_s)*h2o2
			B[1][0] = -1j*r*omk
			B[1][1] = 0
			
			C[0][0] = r*(nub+3*nus)*omk
			
			A[0,:] /= (2*r*omk)
			B[0,:] /= (2*r*omk)
			C[0,:] /= (2*r*omk)
		
			
			M = -(A - 2*C)/(self.dr**2)
			U = -(C/(self.dr**2) + .5*(B-C))/(r*self.dr)
			L = -(C/(self.dr**2) - .5*(B-C))/(r*self.dr)
		
		
			
			
			
			
			
			
			
		
		
		return M,U,L,K