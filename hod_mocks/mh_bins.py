# Find an adeuquate binning in halo mass, 
# such that there are about 0.00025 gal/volume
# or 250000 galaxies in the 1000Mpc/h box.
import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg') 
from matplotlib import pyplot as plt

nmin = 50 #250000

# Directory with the OuterRim simulation haloes
halodir = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/'
istep = 266

outdir = halodir+'ascii/OuterRim_STEP'+str(istep)+'_z0.865/subvols27/'
filename = outdir+'OuterRim_STEP'+str(istep)+'_fofproperties_012.txt'
if (not os.path.exists(filename)):
    print('ERROR: file missing ',filename) ; sys.exit()

# Initial histogram
mmin = 9. ; mmax = 18. ; dmi = 0.05
mbins = np.arange(mmin,mmax,dmi)

mlow, mhigh, mmid = [np.empty((0,1), float) for i in range(3)]
nhini = 0

# Find number of galaxies per initial bin
ff = open(filename) ; iline = 0
for line in ff:
    mm = float(line.split()[6])
    if (mm>0.):
        H, bins_edges = np.histogram(mm,bins=np.append(mbins,mmax))
        nhini = nhini + H

        # Testing -------------
        iline += 1
        if iline>100:
            break
        #-------------------------

mlow = np.append(mlow, mmin)

ii = -1
for jj in range(len(mbins)):
    nn = 0  
    while nn<nmin:
        ii += 1
        nn = nn + nhini[ii]
    print ii,nhini[ii],nn,nmin
    mhigh = np.append(mhigh, mbins[ii])
    mlow = np.append(mlow, mbins[ii])

    if jj<
