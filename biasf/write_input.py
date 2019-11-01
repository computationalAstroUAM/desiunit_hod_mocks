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
	f.write('use_tree= 0 \n')
	f.write('max_tree_order= 6 \n')
	f.write('max_tree_nparts= 100 \n')
	f.write('use_pm= 0 \n')
	f.write('n_grid_side= 256 \n')
	f.close()
        return

def writeScriptRun(runfile, paramfile, outdir):
	f=open(runfile,'w')
	f.write('#!/bin/bash \n')
	f.write('# \n')
	f.write('#SBATCH --nodes=1 \n')
	f.write('#SBATCH --ntasks=1 \n')
	f.write('#SBATCH --time=0-9:00:00  \n')
	f.write('#SBATCH -p sciama4-12.q \n')
	f.write('#SBATCH --job-name=cute  \n')
	f.write('#SBATCH -o '+outdir+'.%j.log \n')
	f.write('# \n')
	f.write('cute=/mnt/lustre/eboss/CUTE/CUTE_box/CUTE_box \n')
	f.write('path='+paramfile+' \n')
	f.write(' \n')
	f.write('$cute $path \n')
	f.write(' \n')
	f.close()
        return

Testing = False

istep = 266 ; iz = 0.865
#241 1.055
#253 0.959
#266 0.865
#279 0.779
#300 0.656

# OuterRim
lbox = 1000.
pathsubvol = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/ascii/OuterRim_STEP'+str(istep)+'_z'+str(iz)+'/subvols27/' 
nameroot = pathsubvol +'OuterRim_STEP'+str(istep)+'_fofproperties_'

halodir = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/'
mfile = halodir+'hmfs/hmf_'+str(istep)+'.txt' 
mlow, mhigh = np.loadtxt(mfile, usecols= (0, 1), unpack=True)
if Testing:
        mlow = np.array([10.58, 10.7]) 

# Output paths
logpath = '/mnt/lustre/gonzalev/Junk/'
pathtemp = logpath+'xi_mass_bin/' ; print('Output: {}'.format(pathtemp))
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
############################## 

namebox = xbox+ybox+zbox
massfileroot = pathtemp + space + '_' + 'full' + '_' + namebox + '_' + 'xyz' + '_'

# Create the run*.sh an parameter files for CUTE
for i, imlow in enumerate(mlow):
        outfile = ocat+space+'_'+'full'+'_'+namebox+'_'+str(imlow)+'.txt'
        massfile = massfileroot + str(imlow) + '.txt'
	if (os.path.exists(massfile)):
		os.remove(massfile)

	paramfile = pathtemp+'paramfile_full_'+namebox+'_'+str(imlow)+'.txt'
	if (os.path.exists(paramfile)):
		os.remove(paramfile)

	runfile = pathtemp+'runfile_full_'+namebox+'_'+str(imlow)+'.sh'
	if (os.path.exists(runfile)):
		os.remove(runfile)

        # Write input parameter file for CUTE
	writeScriptParams(paramfile,massfile,outfile, str(lbox))
        # Write run script for Sciama
	writeScriptRun(runfile, paramfile, logroot)

# Path to file with halo information
boxfile = nameroot + namebox + '.txt'
if (not os.path.exists(boxfile)):
	print ('STOP: This file does not exist {}'.format(boxfile))
        sys.exit()
else:
	print ('Reading: {}'.format(boxfile))

# Initialize temporary arrays and indexes
massh, xx, yy, zz = [np.array([]) for i in range(4)]
ival = 0
val = 50000 

# Get number of lines for the halo file
num_lines = sum(1 for line in open(boxfile)) - 1

# Read file line by line
ff = open(boxfile, 'r')
for iline, line in enumerate(ff):
	massh = np.append(massh, float(line.split()[6]))
	xx = np.append(xx, float(line.split()[0]))
	yy = np.append(yy, float(line.split()[1]))
	zz = np.append(zz, float(line.split()[2]))	

        ival += 1
        # Write files in manageble batches in terms of memory
        if ((ival > val) or (iline == num_lines) ):
		for i, imlow in enumerate(mlow):
                        massfile = massfileroot + str(imlow) + '.txt'
			imhigh = mhigh[i] 

			ind = np.where((massh>=imlow) & (massh<imhigh))
			if (np.shape(ind)[1]>0):
				tofile = zip(xx[ind], yy[ind], zz[ind], massh[ind])
                                if Testing:
                                        print('{} to {}'.format(np.shape(ind)[1],massfile))
				if os.path.exists(massfile):
					with open(massfile, 'a') as outf:
						np.savetxt(outf,tofile,fmt=('%10.5f %10.5f %10.5f %10.5f'))
                                                outf.closed
				else:
					with open(massfile, 'w') as outf:
						np.savetxt(outf,tofile,fmt=('%10.5f %10.5f %10.5f %10.5f'))
                                                outf.closed

		# Reset arrays				
		ival = 0
		massh, xx, yy, zz = [np.array([]) for i in range(4)]
		print ('{} lines read'.format(iline))

        #Testing------------
        if (Testing and iline>3*val):
                print ('Testing finished')
                break
        #-------------------

ff.close()

print ('Program finished')
