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

ifile = 0
for ix in range(ncell):
    for iy in range(ncell):
        for iz in range(ncell):
            filename = root+str(ix)+str(iy)+str(iz)+'.txt'
            if (not os.path.exists(filename)):
                print('ERROR: file missing ',filename)
            else: 
                # Count the lines in the file
                ff = open(filename) ; num = 0
                for line in ff:
                    num += 1
                ff.close()
                print(num,' lines in ',filename)

                # Downsample to plot
                sampling_rate = 0.001
                if num>10000000:
                    sampling_rate = 0.0001

                xp, yp, zp = [np.empty((0,1), float) for i in range(3)]
                ff = open(filename) 
                for line in ff:
                    xx = float(line.split()[0])
                    yy = float(line.split()[1])
                    zz = float(line.split()[2])

                    ran = np.random.random_sample()
                    if (ran < sampling_rate):
                        xp = np.append(xp, xx)
                        yp = np.append(yp, yy)
                        zp = np.append(zp, zz)
                ff.close()                                

                # Plot
                print('Plot ',len(xp),' points')
                fig = plt.figure() 
                ax = fig.add_subplot(111,projection='3d')
                xtit ='x (Mpc/h)' ; ytit ='y (Mpc/h)'; ztit ='z (Mpc/h)'

                ax.set_xlabel(xtit) ; ax.set_ylabel(ytit) ; ax.set_zlabel(ztit)
                ax.scatter(xs=xp,ys=yp,zs=zp)
                
                plotfile = outdir+'plots/OuterRim_STEP'+str(istep)+\
                    '_'+str(ix)+str(iy)+str(iz)+'.pdf'
                fig.savefig(plotfile)
                print 'Test plot: ',plotfile

                # Testing -------------
                ifile += 1
                if ifile>2:
                    break
                #-------------------------
