import sys,os
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

# Read the HMF mass bins
halodir = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/'
mfile = halodir + 'hmf.txt'
mlow, mhigh, mhf, mhmed = np.loadtxt(mfile, usecols= (0, 1, 2, 3), unpack=True)

# Set up the subhistograms:
mmin = np.min(mlow) ; mmax = np.max(mhigh)
dm= 0.01 ; m_arr =  np.arange(mmin,mmax,dm)
mhist = m_arr + dm*0.5

histc = np.zeros(shape=(len(mlow),len(mhist)))
hists = np.zeros(shape=(len(mlow),len(mhist)))

# Mock
#(1) X 	(2) Y 	(3) Z        	position in Mpc/h
#(4) Vx	(5) Vz	(6) Vz 	 	physical velocitie in km/h
#(7) M_h				mass of the host halo in M_sun/h
#(8) Nsat  if the galaxy is satellite: -1; if the galaxy is central: number of satellites associated to the same halo; 
volume = 1000.**3.
mock = '/users/savila/ELGs_eBOSS/HOD_gauss/output/galaxies_1000Mpc_V0.2_mu11.84_Ac0.023_As0.0034_mock000.dat'
mh1, nsat = np.loadtxt(mock,usecols=(6,7), unpack=True)
mh = np.log10(mh1)

for i, imlow in enumerate(mlow):
    imhigh = mhigh[i]

    ind = np.where((mh>=imlow) & (mh<imhigh) & (nsat<0.))
    H, bin_edges = np.histogram(mh[ind],bins=np.append(m_arr,mmax))
    hists[i,:] = H/volume/dm

    ind = np.where((mh>=imlow) & (mh<imhigh) & (nsat>=0.))
    H, bin_edges = np.histogram(mh[ind],bins=np.append(m_arr,mmax))
    histc[i,:] = H/volume/dm

plotname = halodir+'plots/mock_test.png'
fig = plt.figure(figsize = (8., 9.))
xtit = '${\\rm logM_{h}}$' ; ytit = '${\\rm logN}$'
plt.xlabel(xtit) ; plt.ylabel(ytit)
plt.ylim([-6, 0])

plt.plot(mhmed, np.log10(mhf), 'k',label='HMF')

for i, imlow in enumerate(mlow):
    yc = histc[i,:] ; ys = hists[i,:]

    if (i == 1):
        plt.plot(mhist, np.log10(yc), 'g--', label='Centrals')
        plt.plot(mhist, np.log10(ys), 'b:', label='Satellites')
        plt.plot(mhist, np.log10(yc+ys), 'r-', label='All')
    else:
        plt.plot(mhist, np.log10(yc), 'g--')
        plt.plot(mhist, np.log10(ys), 'b:')
        plt.plot(mhist, np.log10(yc+ys), 'r-')

leg = plt.legend(loc=1)
leg.draw_frame(False)
plt.show()
fig.savefig(plotname)
print ('Ouput: {}'.format(plotname))
    
                 
