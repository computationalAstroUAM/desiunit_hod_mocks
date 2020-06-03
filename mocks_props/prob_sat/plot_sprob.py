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

typemock=['NFW'] #,'part']
lsty=['-','-'] #lsty=['-',':']

xboxs = ['0','1','2']
yboxs = ['0','1','2']
zboxs = ['0','1','2']

if (Testing):
	typemock = ['NFW1']
	xboxs = ['0'] ; yboxs = ['1'] ; zboxs =['2']

# Bins in distance (Mpc/h)
nmin = 0.5 ; nmax = 50. ; dn = 1
nbins = np.arange(nmin,nmax,dn)
nhist = nbins +dn*0.5

# Figure 
fig = plt.figure(figsize=(7,4))
xtit = "Number of Satellites"
ytit = "Normalised count"
xmin = 0.9 ; xmax = 5.1

ax = fig.add_subplot(111)
ax.set_xlim(xmin,xmax) 
ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
ind = np.where(nhist<=xmax)
ax.set_xticks(list(nhist[ind]),minor=False)
ax.tick_params(axis='x', which='minor', bottom=False, top=False)

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
maxn = -999.
for itm, tmock in enumerate(typemock):
	# Path to nsat
	nsatpath = '/mnt/lustre/eboss/OuterRim/mocks/nsat_'

	with open(tmock+'_mocks.txt', 'r') as ff:
		mocks = [line.strip() for line in ff]

	if ('1' in tmock):
		tmock = tmock.split('1')[0] 

	for im, mock1 in enumerate(mocks):
		mock2 = mock1.split('/')[-1]
		mock = mock2.replace('galaxies','nsat')
		print('Mock: {}'.format(mock))

		beta1 = mock.split('beta')[1] 
		beta = 'Neg-bin $\\beta$='+beta1.split('_')[0] 
		if (beta1.split('_')[0] == '0.000'):
			beta = 'Poisson'
			lstyle = lsty[itm]
		elif (beta1.split('_')[0] == '-2.000'):
			beta = 'Next integer'
			lstyle = lsty[itm]
		else:
			lstyle = '-' #'--'

		pnsat = np.zeros(shape=(len(nhist)))
		for xbox in xboxs:
			for ybox in yboxs:
				for zbox in zboxs:
					ibox = xbox+ybox+zbox ; print('Box: {}'.format(ibox))

					# Change the mock names to the box we are working with
					imock = mock.replace('mock000','mock'+ibox)
					mockfile = nsatpath+tmock+'/'+imock
					file_fine = check_file(mockfile) 
					if (not file_fine):
						print('Moving on, file not found {}'.format(mockfile))
						break
					#print('Mockfile: {}'.format(mockfile))

					# Read the catalogue with number of satellites
					# tag 0, nsat 1
					lnsat = []
					with open(mockfile, 'rb') as ff:
						for line in ff:
							lnsat.append(int(line.strip().split()[1]))

					nsat = np.asarray(lnsat,dtype=int) 
					lnsat = []
					if (np.max(nsat) > maxn): maxn=np.max(nsat)

					# Histogram
					ind = np.where(nsat>0)
					H, bin_edges = np.histogram(nsat[ind], bins=np.append(nbins,nmax))
					pnsat = pnsat + H
					nsat = []

		print('    Nmax={}'.format(maxn))
		intpn = np.sum(pnsat)*dn
		pn = pnsat/intpn
		# Lines plot
		#ax.plot(nhist,pn,label=tmock+', '+beta,
		#		linestyle=lsty[itm],color=cols[im])
		# Step plot
		tmp = np.insert(pn,0,pn[0])
		yy = np.asarray(tmp,dtype=float)
		tmp = np.insert(nhist,0,nhist[0]-1)
		xx = np.asarray(tmp,dtype=float) + 0.5
		#ax.step(xx,yy,label=tmock+', '+beta,
		#		linestyle=lstyle,color=cols[im])
		ax.step(xx,yy,label=beta,
				linestyle=lstyle,color=cols[im])

# Legend
leg = ax.legend(loc=1)
#leg.draw_frame(False)

# Save figure
plotfile = '/mnt/lustre/eboss/OuterRim/mocks/plots/prob_sat_'+str(istep)+'.png'
#plotfile = '/mnt/lustre/eboss/OuterRim/mocks/plots/prob_sat_bestfit_'+str(istep)+'.png'
fig.savefig(plotfile,dpi=300)
print('Output: ',plotfile)
