import os

def grab_from_makefile(makefile_name):
  
  with open(makefile_name,"r") as f:
    for line in f.readlines():
      if 'SOURCES=' in line:
	files = line.split('=')[-1]
	break
  
  files=files.split()
  return files

def allocation_check(files,check_all=False):
  # This function will try to match every malloc and calloc call with a corresponding free call
 
  allocated_vars=[]
  freed_vars=[]
  allocated_files=[]
  freed_files=[]
  if type(files) == str:
    if len(files) == 0:
      files =[f for f in os.listdir('./') if f[-2:]=='.c']
    else:
      if files[-1]=='/':
        files = [files + f for f in os.listdir(files) if f[-2:]=='.c']
      else:
	files = [files]
  else:
    if type(files) != list:
      print 'Wrong type for input argument'
      return
  
 
    
    
    
  for file_name in files:
    with open(file_name,"r") as f:
      for line in f.readlines():
	if '//' not in line and '/*' not in line and '*/' not in line:
	  if 'malloc(' in line or 'calloc(' in line:
	    commands = line.split(';')
	    for temp in commands:
	      if 'malloc(' in temp or 'calloc(' in temp:
		temp = temp.split('=')[0]
		if '*' in temp:
		  temp = temp.split('*')[-1]
		temp = temp.strip('\t')
		temp = temp.strip('\n')
		temp = temp.strip()  
		allocated_vars.append(temp)
		allocated_files.append((temp,file_name,line))
	  if 'free(' in line and '_free(' not in line:
	    commands = line.split(';')
	    for temp in commands:
	      if 'free(' in temp and '_free(' not in temp:
		temp = (temp.split('free(')[-1]).split(')')[0]
		temp = temp.strip('\t')
		temp = temp.strip('\n')
		temp = temp.strip()
		freed_vars.append(temp)  
		freed_files.append((temp,file_name,line))
    
    
  print 'Allocated Variables:'
  for x in allocated_files:
    print x


  print '\n\nFreed Variables'
  for x in freed_files:
    print x
  

  print '\n\n'

    
  if len(allocated_vars) > len(freed_vars):
    print '%d More allocated variables than freed variables!' % (len(allocated_vars)-len(freed_vars))
      
    for fx in freed_vars:
      if fx in allocated_vars:
	allocated_vars.remove(fx)
      
    print '\nExtra Allocated Variables are:'
    print allocated_vars
      
      
    
  elif len(allocated_vars) < len(freed_vars):
    print '%d More freed variables than allocated variables!' % (len(freed_vars)-len(allocated_vars))
    for ax in allocated_vars:
      if ax in freed_vars:
	freed_vars.remove(ax)
	
    print '\nExtra Freed Variables are:'
    print freed_vars
    
    
  else:
    print 'Same Number for Allocated anad Freed Variables'
      
    allocated_vars = sort(allocated_vars)
    freed_vars = sort(freed_vars)
    for i in range(len(freed_vars)):
      if allocated_vars[i] != freed_vars[i]:
	
	print '\nVariable Mismatch for:'
	print allocated_vars[i], freed_vars[i]