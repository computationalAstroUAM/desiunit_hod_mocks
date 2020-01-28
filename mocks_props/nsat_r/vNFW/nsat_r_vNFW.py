import sys,os
import numpy as np
from iotools import check_file
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec
from distinct_colours import get_distinct
import mpl_style
plt.style.use(mpl_style.style1)

Testing = False

istep = 266

############################
#
# Read K from input params
#
narg = len(sys.argv)
if(narg == 2):
	kmock = sys.argv[1]
else:
	sys.exit('Argument to be passed: K')  
############################

typemock=['NFW','part']
lsty=['-',':']

xboxs = ['0','1','2']
yboxs = ['0','1','2']
zboxs = ['0','1','2']

if (Testing):
	xboxs = ['0'] ; yboxs = ['1'] ; zboxs =['2']

# Bins in distance (Mpc/h)
rmin = 0.01 ; rmax = 5. ; dr = 0.1
rbins = np.arange(rmin,rmax,dr)
rhist = rbins +dr*0.5

# Figure 
fig = plt.figure()
xtit = "$r\,h^{-1}{\\rm Mpc}$"
ytit = "${\\rm log}_{10}(N_{\\rm sat}/N_{\\rm sat, max})$"
xmin = 0. ; xmax = rmax
ymin = 0. ; ymax = 5.

ax = fig.add_subplot(111)
ax.set_xlim(xmin,xmax) #; ax.set_ylim(ymin,ymax)
ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)

# Count mocks to set colour
colmax = 0
for tmock in typemock:
	count = 0
	with open(tmock+'_mocks.txt', 'r') as ff:
		for line in ff: count += 1
	if (count>colmax) : colmax = count
cols = get_distinct(colmax)

# Loop over mocks and boxes
minr = 999. ; maxr = -999.
for itm, tmock in enumerate(typemock):
	# Path to mocks
	mockdir = '/mnt/lustre/savila/HOD_'+tmock+'/output_V1/' 

	# Initialise the matrix to store the number of sat at different distances
	with open(tmock+'_mocks.txt', 'r') as ff:
		mocks = [line.strip() for line in ff]

	for im, mock in enumerate(mocks):
		tmp = mock.replace('_K1.00','_K'+kmock) ; mock = tmp

		if (tmock == 'NFW'):
			vv1 = mock.split('Mpc_')[1]
			vv = vv1.split('more')[0]
			inleg = tmock+', '+vv
			if (im==0): inleg = inleg+', K'+kmock
		else:
			inleg = tmock

		nsatr = np.zeros(shape=(len(rbins)))
		for xbox in xboxs:
			for ybox in yboxs:
				for zbox in zboxs:
					ibox = xbox+ybox+zbox #; print('Box: {}'.format(ibox))

					# Change the mock names to the box we are working with
					imock = mock.replace('mock000','mock'+ibox)
					mockfile = mockdir+imock
					check_file(mockfile) ; print('Mockfile: {}'.format(mockfile))

					# Read mock catalogue
					# x 0, y 1,z 2(Mpc/h), vx 3, vy 4, vz 5(comoving peculiar in km/s), 
					# lmass 6(np.log10(fof_halo_mass)), cen 7(-1=sat, 0=cen), 
					# dvx 8, dvy 9, dvz 10, dx 11, dy 12, dz 13, tag 14(fof_halo_tag)  
					ldx = [] ; ldy = [] ; ldz = []
					with open(mockfile, 'rb') as ff:
						for line in ff:
							cen = int(line.strip().split()[7])
							if (cen<0):
								ldx.append(float(line.strip().split()[11]))
								ldy.append(float(line.strip().split()[12]))
								ldz.append(float(line.strip().split()[13]))
					dx = np.asarray(ldx,dtype=float)
					dy = np.asarray(ldy,dtype=float)
					dz = np.asarray(ldz,dtype=float)
					ldx = [] ; ldy = [] ; ldz = []

					# Get r in kpc/h for satellite galaxies
					rsat = np.sqrt(dx*dx + dy*dy + dz*dz)
					print(np.min(rsat),np.max(rsat))
					if(np.min(rsat)<minr) : minr = np.min(rsat)
					if(np.max(rsat)>maxr) : maxr = np.max(rsat)

					# Histogram
					H, bin_edges = np.histogram(rsat, bins=np.append(rbins,rmax))
					nsatr = nsatr + H
					#print('    Nmin={:.2f}, Nmax={:.2f}'.format(np.min(nsatr),np.max(nsatr)))

		ymax = np.max(nsatr)
		ind = np.where(nsatr>0)
		ax.plot(rhist[ind],np.log10(nsatr[ind]/ymax),label=inleg,
				linestyle=lsty[itm],color=cols[im])

print('    Rmin={:.2f}, Rmax={:.2f} Mpc'.format(minr,maxr))

# Legend
leg = ax.legend(loc=1)
leg.draw_frame(False)

# Save figure
plotfile = '/mnt/lustre/eboss/OuterRim/mocks/plots/nsat_r_K'+kmock+'_'+str(istep)+'.png'
fig.savefig(plotfile)
print 'Output: ',plotfile
