#!/usr/bin/env/python

from sys import argv

if len(argv)>1:	algfile = str(argv[1]) 
else:	algfile = 'explicit' 

with open("Makefile.in","r") as f:
	with open("Makefile","w") as g:
		for line in f.readlines():
			if line.split('=')[0] == 'SOURCES':
				files = line.split('=')[-1].split()
				if algfile == '--implicit':
					files.append('implicit.c')
				else:
					files.append('algo.c')
				
				files = ' '.join(files)
				files += '\n'
				line = '='.join([line.split('=')[0],files])
			g.write(line)
			
			
			