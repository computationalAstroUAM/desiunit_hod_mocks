# Find an adeuquate binning in halo mass, 
# such that there are about 0.00025 gal/volume
# or 250000 galaxies in the 1000Mpc/h box.
import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg') 
from matplotlib import pyplot as plt

nmin = 250000

# Directory with the OuterRim simulation haloes
halodir = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/'
istep = 266

outdir = halodir+'ascii/OuterRim_STEP'+str(istep)+'_z0.865/subvols27/'
filename = outdir+'OuterRim_STEP'+str(istep)+'_fofproperties_012.txt'
if (not os.path.exists(filename)):
    print('ERROR: file missing ',filename) ; sys.exit()

# Initial histogram
mmin = 8. ; mmax = 18. ; dmi = 0.05
mbins = np.arange(mmin,mmax,dmi)

mlow, mhigh, mmid, nbin = [np.empty((0,1), float) for i in range(4)]
nhini = 0

# Find number of galaxies per initial bin
ff = open(filename) ; iline = 0
for line in ff:
    mm = float(line.split()[6])
    if (mm>0.):
        H, bins_edges = np.histogram(mm,bins=np.append(mbins,mmax))
        nhini = nhini + H

        # Testing -------------
        #iline += 1
        #if iline>100:
        #    break
        #-------------------------

ind = np.where(nhini>0)
nhno0 = nhini[ind] ; mbno0 = mbins[ind]

mhigh = np.append(mhigh, mmax) 
nn = 0 
for ii in xrange(len(nhno0)-1, 0, -1):
    nn = nn + nhno0[ii]
    if (nn>=nmin):
        mlow = np.append(mlow, mbno0[ii])
        mhigh = np.append(mhigh, mbno0[ii]) 
        nbin = np.append(nbin, nn)
        nn = 0
mlow = np.append(mlow, mmin)
nbin = np.append(nbin, nn)

print ('Mlow=',mlow)
print ('Mhigh=',mhigh)
print ('nbin=',nbin)

# 50 bins result
#('Mlow=', array([12.85, 12.6 , 12.45, 12.35, 12.25, 12.15, 12.05, 12.  , 11.95,
#       11.9 , 11.85, 11.8 , 11.75, 11.7 , 11.65, 11.6 , 11.55, 11.5 ,
#       11.45, 11.4 , 11.35, 11.3 , 11.25, 11.2 , 11.15, 11.1 , 11.05,
#       11.  , 10.95, 10.9 , 10.85, 10.8 , 10.75, 10.7 , 10.65, 10.6 ,
#        8.  ]))
#('Mhigh=', array([18.  , 12.85, 12.6 , 12.45, 12.35, 12.25, 12.15, 12.05, 12.  ,
#       11.95, 11.9 , 11.85, 11.8 , 11.75, 11.7 , 11.65, 11.6 , 11.55,
#       11.5 , 11.45, 11.4 , 11.35, 11.3 , 11.25, 11.2 , 11.15, 11.1 ,
#       11.05, 11.  , 10.95, 10.9 , 10.85, 10.8 , 10.75, 10.7 , 10.65,
#       10.6 ]))
#('nbin=', array([ 269615.,  286156.,  279933.,  253613.,  319605.,  406039.,
#        511726.,  303355.,  341839.,  386143.,  421567.,  482066.,
#        533405.,  597654.,  682057.,  742721.,  819859.,  934349.,
#       1063111., 1195514., 1228484., 1536951., 1508682., 1880218.,
#       1923083., 2117475., 2646660., 2466211., 3039567., 3847746.,
#       3190423., 3881189., 4844115., 4510534., 5578230., 7114339.,
#             0.]))
