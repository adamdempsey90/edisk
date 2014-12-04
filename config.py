#!/usr/bin/env/python

from sys import argv


with open('inputs/params.opt','r') as f:
	temp = [x.split('+') for x in f.readlines()]
	defs=[]
	for x in temp:
		if x[0] == '' and '#' not in x:
			defs.append(x[-1].split('\n')[0])
	with open('src/defines.h','w') as g:
		if x!= []:
			for x in defs:
				g.write('#define ' + x + '\n')
		else:
			g.write('\n')
print defs
if 'IMPLICIT' in defs:
	algfile = 'implicit.c'
else:
	algfile = 'algo.c rkf.c rk45.c'


if 'INDIRECT' in defs:
	algfile += ' star.c'
	
if 'COMPANION' in defs:
	algfile += ' companion.c'	
with open("Makefile.in","r") as f:
	with open("Makefile","w") as g:
		for line in f.readlines():
			if line.split('=')[0] == 'SOURCES':
				files = line.split('=')[-1].split()
				files.append(algfile)
				files = ' '.join(files)
				files += '\n'
				line = '='.join([line.split('=')[0],files])
			g.write(line)
			
			