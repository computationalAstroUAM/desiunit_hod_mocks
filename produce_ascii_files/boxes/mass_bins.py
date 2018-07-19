# Find an adeuquate binning in halo mass, 
# such that there are about 0.00025 gal/volume
# or 250000 galaxies in the 1000Mpc/h box.
import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg') 
from matplotlib import pyplot as plt

def histedges_equalN(x, nbin): 
    npt = len(x)
    return np.interp(np.linspace(0, npt, nbin + 1),
                     np.arange(npt),
                     np.sort(x))

val = 300000

# Directory with the OuterRim simulation haloes
halodir = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/'
istep = 266

outdir = halodir+'ascii/OuterRim_STEP'+str(istep)+'_z0.865/subvols27/'
root = outdir+'OuterRim_STEP'+str(istep)+'_fofproperties_'

filename = root+'012.txt'
if (not os.path.exists(filename)):
    print('ERROR: file missing ',filename)
else: 
    # Count the lines in the file
    ff = open(filename) ; num = 0
    for line in ff:
        num += 1
    ff.close()
    nbin = int(num/val) 
    if (nbin>50):
        nbin = 50
    nbin = 3
    print('Number of bins set to ',nbin)

    lm = np.empty((0,1), float)
    ff = open(filename) ; iline = 0
    for line in ff:
        mm = float(line.split()[6])
        if (mm>0.):
            lm = np.append(lm, mm)
                               
        # Testing -------------
        #iline += 1
        #if iline>10:
        #    break
        #-------------------------
    ff.close() 

print ('File read')
edges = histedges_equalN(lm, nbin)
print('Edges=',edges)

# Plot
fig = plt.figure() 
ax = fig.add_subplot(111)
xtit ='Mh (Msun/h)' ; ytit ='Number'
ax.set_xlabel(xtit) ; ax.set_ylabel(ytit) 
ax.hist(lm,edges)

plotfile = outdir+'plots/mh_bins_test.pdf'
fig.savefig(plotfile)
print 'Test plot: ',plotfile

