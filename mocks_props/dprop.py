import sys,os.path
import numpy as np

space = 'rspace' #'zspace'
step = 266
redshift = 0.865

# Path to boxes
haloroot = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/ascii/OuterRim_STEP'+\
		  str(step)+'_z'+str(redshift)+'/subvols27/OuterRim_STEP'+\
		  str(step)+'_fofproperties_'
haloend = '.txt'

# Path to mocks
mockdir = '/mnt/lustre/savila/HOD_NFW/output_V1/'

#######################################
#
# Read the box from the input parameters
# 
narg = len(sys.argv)
if(narg == 4):
	xbox = sys.argv[1]
	ybox = sys.argv[2]
	zbox = sys.argv[3]
else:
	sys.exit('Arguments to be passed: xbox, ybox, zbox (0-2)')

#######################################
ibox = xbox+ybox+zbox
halof = haloroot+ibox+haloend
# Check that the file exists
if (not os.path.isfile(halof)): 
	print('STOP: no input file {}'.format(halof))

# Read the mock catalogue file names
ff = open('mocks')
mockfiles = [line.rstrip('\n') for line in ff.readlines() if line.strip()]
ff.close()

# Change the mocks names to the box we are working with
mockfiles = [m.replace('mock000','mock'+ibox) for m in mockfiles]

# 
