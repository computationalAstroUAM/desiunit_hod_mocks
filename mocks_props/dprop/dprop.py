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
print(haloend)
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
ibox = xbox+ybox+zbox ; print('Box: {}'.format(ibox))
halof = haloroot+ibox+haloend
check_file(halof) # Check that the file exists

# Path to mocks
mockdir = '/mnt/lustre/savila/HOD_'+typemock+'/output_V1/'
# Change the mock names to the box we are working with
imock = mock.replace('mock000','mock'+ibox) 
mockfile = mockdir+imock 
check_file(mockfile) ; print('Mockfile: {}'.format(mockfile))

# Read the mock file to memory
if (typemock == 'NFW'): 
	# x 0, y 1,z 2 (Mpc/h), vx 3, vy 4, vz 5 (comoving peculiar in km/s), 
	# lmass 6 (np.log10(fof_halo_mass)), cen 7 (-1=sat, 0=cen), tag 8 (fof_halo_tag) 
	with open(mockfile, 'rb') as ff:
		data = np.array([line.strip().split() for line in ff],float)
elif (typemock == 'part'): 
	#x 0, y 1,z 2(Mpc/h), vx 3, vy 4, vz 5(comoving peculiar in km/s), lmass 6(np.log10(fof_halo_mass)), 
	# cen 7(-1=sat, 0=cen), dvx 8, dvy 9, dvz 10, dx 11, dy 12, dz 13, tag 14(fof_halo_tag) 
	with open(mockfile, 'rb') as ff:
		data = np.array([line.strip().split() for line in ff],float)
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
tmp =  data[:,8] ; gtag = tmp.astype(int)
tmp =  data[:,7] ; cen = tmp.astype(int)

data = [] ; tmp = []

ind = np.where(cen < 0)
if (np.shape(ind)[1]<1):
	print('STOP: no satellites in {}'.format(mockfile)) ; sys.exit()

stag = gtag[ind]
itags = np.squeeze(ind)

# Initialize arrays for differences
dx,dy,dz,dvx,dvy,dvz = (np.zeros(shape=(len(xg))) for i in range(6))
dx[ind] = -999. ; dy[ind] = -999. ; dz[ind] = -999. 
dvx[ind] =-999. ; dvy[ind] =-999. ; dvz[ind]= -999. 

# Read the box line by line and find the populated haloes
count = 0
with open(halof, 'rb') as ff:
	for line in ff:
		# x 0, y 1,z 2 (Mpc/h), vx 3, vy 4, vz 5 (comoving peculiar in km/s), 
		# lmass 6 (np.log10(fof_halo_mass)), tag 7 (fof_halo_tag) 
		htag = int(line.strip().split()[7]) 

		#ind = np.where(gtag == htag)  # Using where
		#if (np.shape(ind)[1]<1): continue

		xh = float(line.strip().split()[0])
		yh = float(line.strip().split()[1])
		zh = float(line.strip().split()[2])
		vxh = float(line.strip().split()[3])
		vyh = float(line.strip().split()[4])
		vzh = float(line.strip().split()[5])

		#dx[ind] = xg[ind] - xh # Using where
		#dy[ind] = yg[ind] - yh
		#dz[ind] = zg[ind] - zh
		#dvx[ind] = vxg[ind] - vxh
		#dvy[ind] = vyg[ind] - vyh
		#dvz[ind] = vzg[ind] - vzh

		irms = [] ; update_tag = 'False'
		for ii, itag in enumerate(itags):
			if(htag == stag[ii]):
				dx[itag] = xg[itag] - xh
				dy[itag] = yg[itag] - yh
				dz[itag] = zg[itag] - zh
				dvx[itag] = vxg[itag] - vxh
				dvy[itag] = vyg[itag] - vyh
				dvz[itag] = vzg[itag] - vzh
				irms.append(ii) # indexes to be removed
				update_tag = 'True'
		
		if (update_tag == 'True' and len(itags)>100):
			for irm in irms:
				tmp = np.delete(stag,irm)
				stag = tmp
				tmp = np.delete(itags,irm)
				itags = tmp

		count +=1
		if((count % 1000000.) == 0): 
			print('{} lines'.format(count)) ; sys.exit()

## Write new files
#outdir = '/mnt/lustre/eboss/OuterRim/mocks/HOD_'+typemock+'/'
#outfile = outdir+imock
#tofile = zip(xg,yg,zg,vxg,vyg,vzg,lmass,cen,
#			 dvx,dvy,dvz,dx,dy,dz,gtag)
#with open(outfile, 'w') as outf:
#	np.savetxt(outf,tofile,fmt=('%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %i %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %i'))
#	outf.closed
#print('Output: {}'.format(outfile))
