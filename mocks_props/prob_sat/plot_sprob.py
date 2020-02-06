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

Testing = True

istep = 266

typemock=['NFW','part']
lsty=['-',':']

xboxs = ['0','1','2']
yboxs = ['0','1','2']
zboxs = ['0','1','2']

if (Testing):
	typemock = ['part']
	xboxs = ['0'] ; yboxs = ['1'] ; zboxs =['2']

# Bins in distance (Mpc/h)
nmin = 0.5 ; nmax = 10.5 ; dn = 1
nbins = np.arange(nmin,nmax,dn)
nhist1 = nbins +dn*0.5
tmp = np.insert(nhist1,0,0)
nhist = np.asarray(tmp,dtype=int)

# Figure 
fig = plt.figure()
xtit = "Number of Satellites"
ytit = "$P(N)$"
xmin = 0. ; xmax = nmax
ymin = 0. ; ymax = 1.

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
maxn = -999.
for itm, tmock in enumerate(typemock):
	# Path to mocks
	mockdir = '/mnt/lustre/savila/HOD_'+tmock+'/output_V1/' 

	# Initialise the matrix to store the number of sat at different distances
	with open(tmock+'_mocks.txt', 'r') as ff:
		mocks = [line.strip() for line in ff]

	for im, mock in enumerate(mocks):
		beta1 = mock.split('beta')[1]
		beta = 'beta='+beta1.split('_')[0]
		if (beta1.split('_')[0] == '0.000'):
			beta = 'Poisson'
		elif (beta1.split('_')[0] == '-2.000'):
			beta = 'Next integer'

		print('Mock: {}'.format(mock))
		pnsat = np.zeros(shape=(len(nbins)))
		nzeros = 0
		for xbox in xboxs:
			for ybox in yboxs:
				for zbox in zboxs:
					ibox = xbox+ybox+zbox #; print('Box: {}'.format(ibox))

					# Change the mock names to the box we are working with
					imock = mock.replace('mock000','mock'+ibox)
					mockfile = mockdir+imock
					check_file(mockfile) #; print('Mockfile: {}'.format(mockfile))

					# Read mock catalogue
					# x 0, y 1,z 2(Mpc/h), vx 3, vy 4, vz 5(comoving peculiar in km/s), 
					# lmass 6(np.log10(fof_halo_mass)), cen 7(-1=sat, 0=cen), 
					# dvx 8, dvy 9, dvz 10, dx 11, dy 12, dz 13, tag 14(fof_halo_tag)  
					ltags = [] ; lstag = [] ; lnsat = []
					with open(mockfile, 'rb') as ff:
						for line in ff:
							if (len(line.strip().split()) == 9):
								itag = 8
							elif (len(line.strip().split()) == 15):
								itag = 14
							else:
								print('STOP: unknown file set-up') ; sys.exit()

							thistag = int(line.strip().split()[itag])
							ltags.append(thistag)

							cen = int(line.strip().split()[7])
							if (cen<0):
								lstag.append(thistag)

					# Check that there are more tags than sat. tags
					if (len(ltags)<len(lstag)):
						print('STOP: number tags={}, n sat.={}'.format(
							len(ltags),len(lstag)))
						sys.exit()
					tags = np.unique(np.asarray(ltags,dtype=int))
					nzeros = nzeros + len(ltags)-len(lstag)

					uniqsat = set(lstag)
					count = 0  #Option1
					for usat in uniqsat:
						lnsat.append(lstag.count(usat))
						count += 1
						if((count % 10000)==0): print(count)
					nsat = np.asarray(lnsat,dtype=int)

					ltags = [] ; lstag = [] ; lnsat = []
					
					if (np.max(nsat)>maxn): maxn = np.max(nsat)
					
					# Histogram
					H, bin_edges = np.histogram(nsat, bins=np.append(nbins,nmax))
					pnsat = pnsat + H
					
		tmp = np.insert(pnsat,0,nzeros)
		pnsat = tmp
		intpn = np.sum(pnsat)*dn

		pn = pnsat/intpn
		ax.plot(nhist,pn,label=tmock+', '+beta,
				linestyle=lsty[itm],color=cols[im])

print('    Nmax={}'.format(maxn))

# Legend
leg = ax.legend(loc=1)
leg.draw_frame(False)

# Save figure
plotfile = '/mnt/lustre/eboss/OuterRim/mocks/plots/prob_sat_'+str(istep)+'.png'
fig.savefig(plotfile)
print 'Output: ',plotfile
