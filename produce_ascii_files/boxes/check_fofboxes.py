import sys, os
import numpy as np
import random
import matplotlib
matplotlib.use('Agg') 
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

istep = 241
iz = 1.055

# 266 0.865 Done
# 241 1.055 Done
# 253 0.959 
# 279 0.779  
# 300 0.656  

# Directory with the OuterRim simulation haloes
halodir = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/'

outdir = halodir+'ascii/OuterRim_STEP'+str(istep)+'_z'+str(iz)+'/subvols27/'
root = outdir+'OuterRim_STEP'+str(istep)+'_fofproperties_'

# Directory with plots
plotdir = outdir+'plots/'
if not os.path.exists(plotdir):
    os.makedirs(plotdir)

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
                print('### ERROR: file missing {}'.format(filename))
            else: 
                #ix = 0 ; iy = 2 ; iz = 2
                # Count the lines in the file
                ff = open(filename) ; num = 0
                for line in ff:
                    num += 1
                ff.close()
                print('{} {} lines in  {}'.format(ifile,num,filename))

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
                print('Plot {} points'.format(len(xp)))
                fig = plt.figure() 
                ax = fig.add_subplot(111,projection='3d')
                xtit ='x (Mpc/h)' ; ytit ='y (Mpc/h)'; ztit ='z (Mpc/h)'

                ax.set_xlabel(xtit) ; ax.set_ylabel(ytit) ; ax.set_zlabel(ztit)
                ax.scatter(xs=xp,ys=yp,zs=zp)
                
                plotfile = plotdir+'OuterRim_STEP'+str(istep)+\
                    '_'+str(ix)+str(iy)+str(iz)+'.pdf'
                fig.savefig(plotfile)
                plt.close(fig)
                print ('Test plot: {} '.format(plotfile))

                # Testing -------------
                #ifile += 1
                #if ifile>1:
                #    sys.exit()
                #-------------------------
