#!/usr/bin/env/python
from subprocess import call

input_directory = 'inputs/'
output_directory = 'outputs/'

with open(input_directory + 'params.in','r') as f:
	lines = f.readlines()


beta_vals = [-1.5, -.75, -.5, 0]
salpha_vals = [ .03, .003, .0003,.00003]
balpha_vals = [ -.02, -.002, -.0002,-.00002]


for i in range(len(beta_vals)):
	for j in range(len(salpha_vals)):
		print i,j
		b = beta_vals[i]
		a = salpha_vals[j]
		
		dir_name = output_directory + 'analytic_'+str(b)+'_'+str(a)
			
		call(['mkdir',dir_name])
	
		for k,line in enumerate(lines):
			if line.split('=')[0] == 'alpha shear ':
				lines[k] = 'alpha shear = ' + str(a) + '\n'
				print lines[i]
			if line.split('=')[0] == 'alpha bulk ':
				lines[k] = 'alpha bulk = ' + str(balpha_vals[j]) + '\n'
			if line.split('=')[0] == 'sigma index ':
				lines[k] = 'sigma index = ' + str(b) + '\n'
		with open(input_directory + 'params.in','w') as f:
			for line in lines:
				f.write(line)
		call(['./edisk'])
		call(['cp',output_directory+'output_0.dat',dir_name])
		call(['cp',output_directory+'output_1.dat',dir_name])
		call(['cp',output_directory+'disk.dat',dir_name])
		call(['cp',output_directory+'params.txt',dir_name])
	
