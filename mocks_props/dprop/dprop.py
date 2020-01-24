import sys,os.path
import numpy as np
import stats
from iotools import check_file

space = 'rspace' #'zspace'
step = 266
redshift = 0.865

# Path to boxes
haloroot = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/ascii/OuterRim_STEP'+\
		  str(step)+'_z'+str(redshift)+'/subvols27/OuterRim_STEP'+\
		  str(step)+'_fofproperties_'
haloend = '.txt'

#######################################
#
# Read the box from the input parameters
# 
narg = len(sys.argv)
if(narg == 6):
	xbox = sys.argv[1]
	ybox = sys.argv[2]
	zbox = sys.argv[3]
	typemock = sys.argv[4] ; typemock = typemock.replace('test','')
	mock = sys.argv[5]
else:
	sys.exit('Arguments to be passed: xbox, ybox, zbox (0-2), typemock, mock')

#######################################
ibox = xbox+ybox+zbox
halof = haloroot+ibox+haloend
check_file(halof) # Check that the file exists

# Path to mocks
mockdir = '/mnt/lustre/savila/HOD_'+typemock+'/output_V1/'
# Change the mock names to the box we are working with
imock = mock.replace('mock000','mock'+ibox)
mockfile = mockdir+imock
check_file(mockfile)

# Read the mock file to memory
if (typemock == 'NFW'): 
	# x 0, y 1,z 2 (Mpc/h), vx 3, vy 4, vz 5 (comoving peculiar in km/s), 
	# lmass 6 (np.log10(fof_halo_mass)), cen 7 (-1=sat, 0=cen), tag 8 (fof_halo_tag) 
	with open(mockfile) as ff:
		data = np.array([line.strip().split() for line in ff],float)
elif (typemock == 'part'): 
	#x 0, y 1,z 2(Mpc/h), vx 3, vy 4, vz 5(comoving peculiar in km/s), lmass 6(np.log10(fof_halo_mass)), 
	# cen 7(-1=sat, 0=cen), dvx 8, dvy 9, dvz 10, dx 11, dy 12, dz 13, tag 14(fof_halo_tag) 
	with open(mockfile) as ff:
		data1 = np.array([line.strip().split() for line in ff],float)
		if (np.shape(data)[1] > 9):
			print('Finish, as file {} seems to already contain differences'.format(mockfile))
			sys.exit()
xg =   data[:,0]
yg =   data[:,1]
zg =   data[:,2]
vxg =  data[:,3]
vyg =  data[:,4]
vzg =  data[:,5]
lmass =data[:,6]
tmp =  data[:,7] ; cen = tmp.astype(int)
tmp =  data[:,8] ; gtag = tmp.astype(int)
data = [] ; tmp = []

# Read the box line by line and find the populated haloes
with open(halof) as ff: 
	# x 0, y 1,z 2 (Mpc/h), vx 3, vy 4, vz 5 (comoving peculiar in km/s), 
	# lmass 6 (np.log10(fof_halo_mass)), tag 7 (fof_halo_tag) 
	data = np.array([line.strip().split() for line in ff],float) 

	xh = data[:,0]
	yh = data[:,1]
	zh = data[:,2]
	vxh =data[:,3]
	vyh =data[:,4]
	vzh =data[:,5]
	tmp =data[:,5] ; htag = tmp.astype(int)
	data = [] ; tmp = []

# Initialize arrays for differences
dx,dy,dz,dvx,dvy,dvz = (np.zeros(shape=(len(xg))) for i in range(6))
dx.fill(-999.) ; dy.fill(-999.) ; dz.fill(-999.) 
dvx.fill(-999.) ; dvy.fill(-999.) ; dvz.fill(-999.) 

for ii,itag in enumerate(gtag):
	ind = np.where(htag == itag)
	if(np.shape(ind)[1] != 1): 
		print('WARNING: tag {} from mock {} not match in halo file {}'.format(itag,mockfile,halof))
		continue

	dx[ii] = xg[ii] - xh[ind]
	dy[ii] = yg[ii] - yh[ind]
	dz[ii] = zg[ii] - zh[ind]
	dvx[ii] = vxg[ii] - vxh[ind]
	dvy[ii] = vyg[ii] - vyh[ind]
	dvz[ii] = vzg[ii] - vzh[ind]

xh = [] ; yh = [] ; zh = [] ; vxh = [] ; vxh = [] ; vxh = [] 

# Write new files
outdir = '/mnt/lustre/eboss/OuterRim/mocks/HOD_'+typemock+'/'
outfile = outdir+imock
tofile = zip(xg,yg,zg,vxg,vyg,vzg,lmass,cen,
			 dvx,dvy,dvz,dx,dy,dz,gtag)
with open(outfile, 'w') as outf:
	np.savetxt(outf,tofile,fmt=('%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %i %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %i'))
	outf.closed
print('Output: {}'.format(outfile))
