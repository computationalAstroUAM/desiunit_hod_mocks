import numpy as np
import sys, os
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

snap = '266'
space = 'rspace' #'zspace'

# Path to CUTE output and output plots
ocat = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/ascii/OuterRim_STEP266_z0.865/subvols27/biasf/'
path2plot = ocat+'plots/'

# Read the HMF mass bins
halodir = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/'
mfile = halodir + 'hmfs/hmf_'+snap+'.txt'

mlow, mhigh = np.loadtxt(mfile, usecols= (0, 1), unpack=True)

# Check CUTE output
check_output = False
if check_output:
	for ix in range(3):
		for iy in range(3):
			for iz in range(3):
				namebox = str(ix)+str(iy)+str(iz)

				for i, imlow in enumerate(mlow):
					boxfile = ocat+space+'_full_'+namebox+'_'+str(imlow)+'.txt'
					if (not os.path.exists(boxfile)):
						print ('This output does not exist: {}'.format(boxfile))
						continue
					try: 
						r, xi, err, dd = np.loadtxt(boxfile, unpack=True)
					except:
						print ('Problem reading r,xi,error,dd from: {}'.format(boxfile))


# Plot CUTE outpus
doplots = True
if doplots:
	for i, imlow in enumerate(mlow):
		fig = plt.figure(figsize = (8., 9.))
		xtit = '${\\rm log_{10}(r /h^{-1}\\rm{Mpc}))}$'
		ytit = '${\\rm log_{10}(\\xi(r))}$'
		plt.xlabel(xtit) ; plt.ylabel(ytit)
		plotname = '2PCF_'+str(imlow)+'.png'

		icount = -1
		for ix in range(3):
			for iy in range(3):
				for iz in range(3):
					icount += 1

					namebox = str(ix)+str(iy)+str(iz)
					boxfile = ocat+space+'_full_'+namebox+'_'+str(imlow)+'.txt'
					r, xi, err, dd = np.loadtxt(boxfile, unpack=True)

					ind = np.where(xi>0.)
					if (icount<1):
						leg = str(imlow)+\
						    '$\leq M_{\\rm halo}(h^{-1}\\rm{M_{\odot}})<$'+\
						    str(mhigh[i])
						plt.plot(np.log10(r[ind]), np.log10(xi[ind]), \
								 label=leg, color = 'k')
					else:
						plt.plot(np.log10(r[ind]), np.log10(xi[ind]), \
								 color = 'k')
		leg = plt.legend(loc=1)
		leg.draw_frame(False)
		fig.savefig(path2plot + plotname) ; plt.close()
		print ('Ouput: {}{}'.format(path2plot,plotname))

