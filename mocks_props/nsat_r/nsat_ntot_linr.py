import sys,os
import numpy as np
from iotools import check_file
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec
#from distinct_colours import get_distinct
#import mpl_style
#plt.style.use(mpl_style.style1)

Testing = False

istep = 266

typemock=['NFW','part']
lsty=['-',':']

xboxs = ['0','1','2']
yboxs = ['0','1','2']
zboxs = ['0','1','2']

if (Testing):
	xboxs = ['0'] ; yboxs = ['1'] ; zboxs =['2']

# Bins in distance (Mpc/h)
rmin = 0.03 ; rmax = 5. ; dr = 0.05
rbins = np.arange(rmin,rmax,dr)
rhist = rbins +dr*0.5

# Figure 
fig = plt.figure(figsize=(7,4.2))
xtit = "$r\,[{\\rm Mpc}h^{-1}]$"
ytit = "${\\rm log}_{10}(N_{\\rm sat,dex}/N_{\\rm sat, tot})$"
xmin = rmin ; xmax = rmax
ymin = -5 ; ymax = -0.4

ax = fig.add_subplot(111)
ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax)
ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)

## Count mocks to set colour
#colmax = 0
#for tmock in typemock:
#	count = 0
#	with open(tmock+'_mocks.txt', 'r') as ff:
#		for line in ff: count += 1
#	if (count>colmax) : colmax = count
#cols = get_distinct(colmax)

# Santi's colours                                                                     
cols = ['#bcbd22','blue','red','green'] 

# Loop over mocks and boxes
minr = 999. ; maxr = -999.
for itm, tmock in enumerate(typemock):
	# Initialise the matrix to store the number of sat at different distances
	with open(tmock+'_mocks.txt', 'r') as ff:
		mocks = [line.strip() for line in ff]

	for im, mock in enumerate(mocks):
		avgr = 0. ; ntot = 0.
		kval1 = mock.split('K')[1]
		kval = kval1.split('_')[0]

		nsatr = np.zeros(shape=(len(rbins)))
		for xbox in xboxs:
			for ybox in yboxs:
				for zbox in zboxs:
					ibox = xbox+ybox+zbox #; print('Box: {}'.format(ibox))

					# Change the mock names to the box we are working with
					mockfile = mock.replace('mock000','mock'+ibox)
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

					# Average r
					avgr = avgr + np.sum(rsat)
					ntot = ntot + len(rsat)

					# Limits of r
					if(np.min(rsat)<minr) : minr = np.min(rsat)
					if(np.max(rsat)>maxr) : maxr = np.max(rsat)

					# Histogram
					H, bin_edges = np.histogram(rsat, bins=np.append(rbins,rmax))
					nsatr = nsatr + H
					print('    Nmin={:.2f}, Nmax={:.2f}'.format(np.min(nsatr),np.max(nsatr)))

		ytot = np.sum(nsatr)
		ind = np.where(nsatr>0)
		ax.plot(rhist[ind],np.log10(nsatr[ind]/ytot),label=tmock+', K='+kval,
				linestyle=lsty[itm],color=cols[im]) 
		#ax.plot(np.log10(rhist[ind]),nsatr[ind]/ymax,label=tmock+', K='+kval,
		#		linestyle=lsty[itm],color=cols[im]) #lineal 

		# Arrows/lines to show the median values
		dy = 0.05*(ymax-ymin) 
		avgr = avgr/ntot
		ax.arrow(avgr,ymax-dy,0,dy,color=cols[im],ls=lsty[itm],
				 length_includes_head=True)
				 #head_width=dy*0.1, head_length=dy*0.3)

		ax.set_xscale('log')

print('    Rmin={:.2f}, Rmax={:.2f} Mpc'.format(minr,maxr))

# Legend
leg = ax.legend(loc=0)
#leg.draw_frame(False)

# Save figure
#plotfile = '/mnt/lustre/eboss/OuterRim/mocks/plots/nsat_r_'+str(istep)+'.png'
plotfile = '/mnt/lustre/eboss/OuterRim/mocks/plots/nsat_ntot_linr_'+str(istep)+'.png'
fig.savefig(plotfile,dpi=300)
print('Output: ',plotfile)
