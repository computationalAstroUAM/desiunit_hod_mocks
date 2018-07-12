import sys, os
import numpy as np
import random
import matplotlib
matplotlib.use('Agg') 
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Directory with the OuterRim simulation haloes
halodir = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/'
istep = 266

outdir = halodir+'ascii/OuterRim_STEP'+str(istep)+'_z0.865/subvols27/'
root = outdir+'OuterRim_STEP'+str(istep)+'_fofproperties_'

# Subvolumes
lbox = 3000.
cell = lbox/3. 
ncell = int(lbox/cell) 

# Remove files if they already exist
for ix in range(ncell):
    for iy in range(ncell):
        for iz in range(ncell):
            filename = root+str(ix)+str(iy)+str(iz)+'.txt'
            if (not os.path.exists(filename)):
                print('ERROR: file missing ',filename)
            else: # Try to read
                xc,yc,zc,vx,vy,vz,lmass,tag = \
                    np.loadtxt(filename, unpack = True)

                num = len(xc)
                if (num<1):
                    print('ERROR: empty file ',filename)
                else:
                    # Downsample to plot
                    idx = np.arange(num) ; random.shuffle(idx)
                    val = num/100

                    # Plot
                    fig = plt.figure()
                    ax = fig.add_subplot(111,projection='3d')
                    xtit ='x (Mpc/h)' ; ytit ='y (Mpc/h)'; ztit ='z (Mpc/h)'
                    ax.set_xlabel(xtit) ; ax.set_ylabel(ytit) ; ax.set_zlabel(ztit)
                    ax.scatter(xs=xc[idx[:val]],ys=yc[idx[:val]],zs=zc[idx[:val]])

                    plotfile = outdir+'plots/OuterRim_STEP'+str(istep)+\
                        '_'+str(ix)+str(iy)+str(iz)+'.pdf'
                    fig.savefig(plotfile)
                    print 'Test plot: ',plotfile

