import sys, os, glob, getpass
sys.path.append('/mnt/lustre/eboss/OuterRim/genericio/python/')
import genericio as gio
import numpy as np
import matplotlib ; matplotlib.use('Agg')                                    
from matplotlib import pyplot as plt  

def percentiles(val,data,weights=None):
	if (val <0 or val >1):
		sys.exit('STOP percentiles: 0<val<1')

	if (weights is None):
		ws = np.zeros(shape=(len(data))) ; ws.fill(1.)
	else:
		ws = weights

	data = np.array(data) ; ws = np.array(ws)
	ind_sorted = np.argsort(data)  # Median calculation from wquantiles
	sorted_data = data[ind_sorted] ; sorted_weights = ws[ind_sorted]

	num = np.cumsum(sorted_weights) - 0.5*sorted_weights 
	den = np.sum(sorted_weights) 
	if (den!=0): 
		pn = num/den   
		percentiles = np.interp(val, pn, sorted_data)  
	else:
		sys.exit('STOP percentiles: problem with weights')

	return percentiles

# Read snapshot
narg = len(sys.argv)
if(narg == 2):
	istep = int(sys.argv[1])
	iz = 0
else:
	sys.exit('1 argument to be passed: step(snapshot number)')

# Make a plot?
makeplot = True
withlegend = False
if makeplot:
	from distinct_colours import get_distinct 
	cols = get_distinct(10) 

# Bins
#edges = np.array([10.58, 10.7, 10.8, 10.9, 11., 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12., 12.1, 12.2, 12.3, 12.4, 12.5, 12.7, 13., 13.5, 14., 16.])
edges = np.concatenate( [ np.array( np.arange(10.58,12.4999,0.04)), np.array(np.arange(12.5,13.6,0.05)),  np.array(np.arange(13.6,14.0,0.1)),np.array(np.arange(14.2,14.6,0.2)), np.array([16.]) ])
print("edges={}".format(edges))

dm = edges[1:]-edges[:-1]
mhist = edges[1:]-0.5*dm 
print("mhist={}".format(mhist))

# Bins for the calculation of the median mass in each bin
mhmed = np.zeros_like(mhist)
nbmed = 100
xbmed = np.zeros((len(mhist),nbmed))
for i,ed in enumerate(edges[:-1]):
	xbmed[i,:] = np.linspace(ed, edges[i+1], num=nbmed)
ybmed = np.zeros((len(mhist),nbmed-1))

# Mean mass in each bin
xmean = np.zeros(len(mhist))

###################
# Directory with the OuterRim simulation haloes
halodir = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/'

# OuterRim simulation characteristics (FOF b=0.168 here)
mp  = 1.9E+09 # Msol/h
lbox= 3000.   # Mpc/h

# Get the conversion between the name of the time step and redshift
step = np.genfromtxt(halodir+'step_redshift.txt',usecols=0,dtype=int)
redshift = np.genfromtxt(halodir+'step_redshift.txt',usecols=1)

# Initialize the parameters for the figure ------------------
plt.rcParams['legend.numpoints'] = 1 
plt.rcParams['axes.labelsize'] = 10.0 ; fs = 20 
plt.rcParams['lines.linewidth'] = 2 
fig = plt.figure(figsize=(8.5,9.))

xtit = "${\\rm log}_{10}(\\rm{M/M_{\odot}}h^{-1})$"
ytit = "${\\rm log}_{10}(\Phi/ Mpc^{-3}h^3 {\\rm dlog}_{10}M)$"

xmin = 10. ; xmax = 16.
ymin = -6.5 ; ymax = 0.

jj = 111
ax = fig.add_subplot(jj) 
ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax)                              
ax.set_xlabel(xtit,fontsize = fs) ; ax.set_ylabel(ytit,fontsize = fs)        
#-------------------------------------------------------------                

# Paths to files
zz = redshift[np.where(step == istep)]
print('Processing snapshot at redshift={}'.format(zz))
nroot = halodir+'HaloCatalog/STEP'+str(istep)    

# Initialize to 0 the halo mass functions
ycount, yh = [np.zeros(len(mhist)) for _ in range(2)]

# Loop over each of the sub volumes
infiles = glob.glob(nroot+'/*'+str(istep)+'.fofproperties#*')
for inf,infile in enumerate(infiles):
	print(infile)
	# Print out once the information stored in each file
	if (iz == 0 and inf == 0):
		gio.gio_inspect(infile)
		print('----------------------------')

	# Read the number of particles per halo
	in_count = mp*gio.gio_read(infile, "fof_halo_count")
	ind = np.where(in_count >0.)
	count = np.log10(in_count[ind])
	H, bins_edges = np.histogram(count,bins=edges)
	ycount = ycount + H

	in_count = [] ; ind = [] # Flush arrays

	# FOF mass (Msun/h)
	in_mh = gio.gio_read(infile, "fof_halo_mass")
	ind = np.where(in_mh >0.) ; mh = np.log10(in_mh[ind])
	H, bins_edges = np.histogram(mh,bins=edges)
	yh = yh + H

	in_mh = [] ; ind = [] # Flush arrays

	# Mean
	H, bins_edges = np.histogram(mh,bins=edges,weights=mh)
	xmean = xmean + H

	# Build up histograms for the calculation of the median
	for i,ed in enumerate(edges[:-1]):
		edmed = xbmed[i,:]
		H, bins_edges = np.histogram(mh,bins=edmed)
		ybmed[i,:] = ybmed[i,:] + H

	mh = []  # Flush array
	# Testing----------------------
	#if inf>1:
	#    break
	#-----------------------------

xmean = xmean/yh
ycount = ycount/dm/(lbox**3)
yh = yh/dm/(lbox**3)

# Calculate the median mass of each HMF bin
for i,ed in enumerate(edges[:-1]):
	x = (xbmed[i,:-1]+xbmed[i,1:])/2.
	y = ybmed[i,:]
	mhmed[i] = percentiles(0.5,x,weights=y)

# Write Halo Mass function to a file
tofile = zip(edges[:-1],edges[1:], ycount, mhmed, mhist, dm, xmean)
outfile = halodir+'hmfs/hmf_'+str(istep)+'.txt'
#outfile = 'hmf.txt'
with open(outfile, 'w') as outf:
	outf.write('# log10Mmin_bin log10Mmax_bin Nhalos/dm/vol log10Mmedian_bin (Mass bin mid point) dM log10Mmean_bin \n')
	np.savetxt(outf,tofile,fmt=('%10.5f %10.5f %6.4e %10.5f %10.5f %10.5f %10.5f'))      
outf.closed

# Plot for this redshift
if makeplot:
	ind = np.where(ycount >0.)
	ax.plot(mhist[ind],np.log10(ycount[ind]),\
			color=cols[iz],label='z='+str(zz))

	ind = np.where(yh >0.)
	ax.plot(mhist[ind],np.log10(yh[ind]),\
		color=cols[len(cols)-iz-1],linestyle=':',label=' from mass')

# Directory with outputs (it'll be generated if it doesn't exist)
outdir = halodir+'hmfs/plots/'
if (not os.path.exists(outdir)): os.makedirs(outdir)
print('Output: {}'.format(outfile))

if makeplot:
	if withlegend:
		# Legend
		plt.legend(loc=3,prop={'size':(fs-2)}) 
		# Save figure
		plotfile = outdir+'outerrim_hmf_'+str(istep)+'.pdf'
	else:
		plotfile = outdir+'outerrim_hmf_'+str(istep)+'_nolegend.pdf'

	fig.savefig(plotfile)
	print('Plot: {}'.format(plotfile))
