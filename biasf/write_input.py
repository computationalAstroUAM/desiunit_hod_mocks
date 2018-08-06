import numpy as np
import os, sys
import glob

def writeScriptParams(paramfile, infile, outfile, lbox):
	f=open(paramfile,'w')
	f.write('data_filename= '+infile+'\n')
	f.write('num_lines= all \n')
	f.write('input_format= 0 \n')
	f.write('output_filename= '+outfile+'\n')
	f.write('box_size= '+lbox+'\n')
#	f.write(' \n')
	f.write('use_tree= 0 \n')
	f.write('max_tree_order= 6 \n')
	f.write('max_tree_nparts= 100 \n')
#	f.write(' \n')
	f.write('use_pm= 0 \n')
	f.write('n_grid_side= 256 \n')
#	f.write(' \n')
	f.close()

def writeScriptRun(runfile, paramfile, outdir):
	f=open(runfile,'w')
	f.write('#!/bin/bash \n')
	f.write('#PBS -l walltime=24:00:00 \n')
	logid = outdir+'.o.${PBS_JOBID}'
	f.write('#PBS -o '+logid+' \n')
	f.write('#PBS -e '+outdir+'.e.${PBS_JOBID} \n')
	f.write('#PBS -l nodes=1:ppn=1 \n')
	f.write('#PBS -q sciama1.q \n')
	f.write(' \n')
	f.write('echo '+logid+' >> jobs_id.txt \n')
	f.write('cute=/users/ghthomas/CUTE/CUTE_box/CUTE_box \n')
#	f.write('cute=/users/ghthomas/outerrim_mocks/george_thomas_project/test.py \n')
	f.write('path='+paramfile+' \n')
	f.write(' \n')
	f.write('$cute $path \n')
#	f.write('python $cute \n')
	f.write(' \n')
	f.close()

# OuterRim
lbox = 1000.
pathsubvol = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/ascii/OuterRim_STEP266_z0.865/subvols27/'
nameroot = pathsubvol + 'OuterRim_STEP266_fofproperties_'

halodir = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/'
mfile = halodir + 'hmf.txt'
mlow, mhigh = np.loadtxt(mfile, usecols= (0, 1), unpack=True)
#mlow = np.array([10.58, 10.7]) # For testing

# Output paths
logpath = '/mnt/lustre/gonzalev/Junk/'
pathtemp = logpath+'xi_mass_bin/' ; print('Output:',pathtemp)
logroot = logpath+'CUTErun'

ocat = pathsubvol+'biasf/'

#############################                                      
#                                                                  
#  Input ARGUMENTS                                                 
#                                                                  
narg = len(sys.argv)                                               
if(narg == 5):
    xbox = sys.argv[1]
    ybox = sys.argv[2]
    zbox = sys.argv[3]
    space = sys.argv[4]
else:                                                              
    sys.exit('4 arguments to be passed: ix, iy, iz, space')

############################# 

namebox = xbox+ybox+zbox

# Create the run*.sh an parameter files for CUTE
for i, imlow in enumerate(mlow):
	outfile = ocat+space+'_'+'full'+'_'+namebox+'_'+str(imlow)+'.txt'

	massfile = pathtemp + space + '_' + 'full' + '_' + namebox + '_' + 'xyz' + '_' + str(imlow) + '.txt'
	if (os.path.exists(massfile)):
		os.remove(massfile)

	paramfile = pathtemp+'paramfile_full_'+namebox+'_'+str(imlow)+'.txt'
	if (os.path.exists(paramfile)):
		os.remove(paramfile)

	runfile = pathtemp+'runfile_full_'+namebox+'_'+str(imlow)+'.sh'
	if (os.path.exists(runfile)):
		os.remove(runfile)
		
	writeScriptParams(paramfile,massfile,outfile, str(lbox))
	writeScriptRun(runfile, paramfile, logroot)
			

# Path to file with halo information
boxfile = nameroot + namebox + '.txt'
if (not os.path.exists(boxfile)):
	print ('This file does not exist,', boxfile)
else:
	print ('Reading:', boxfile)

# Read file line by line
ff =open(boxfile, 'r') ; iline =0
for line in ff:
	massh = float(line.split()[6])

	for i, imlow in enumerate(mlow):
		imhigh = mhigh[i] #; print (imlow, massh, imhigh)

		if ((massh>=imlow) & (massh<imhigh)):
			x = (line.split()[0])
			y = (line.split()[1])
			z = (line.split()[2])

			# Creating the input for CUTE
			massfile = pathtemp+space+'_'+'full'+'_'+namebox+'_'+'xyz'+'_'+str(imlow)+'.txt'
			if os.path.exists(massfile):
				with open(massfile, 'a') as outf:
					outf.write(x+' '+y+' '+z+' \n')
					outf.closed
					
			else:
				with open(massfile, 'w') as outf:
					outf.write(x+' '+y+' '+z+' \n')
					outf.closed
					continue

	##Testing------------
        #iline += 1 #; print (iline)
	#if (iline>100):
	#	print ('Testing finished')
	#	break
	##-------------------

ff.close()

print ('Program finished')
