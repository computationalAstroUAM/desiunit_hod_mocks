import sys,os
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec
#import mpl_style
#plt.style.use(mpl_style.style1)
plt.rcParams.update({'font.size': 16})

def chi2(obs,model,err):
	val = 0.
	for i,iobs in enumerate(obs):
		val = val + (iobs-model[i])**2/(err[i]*err[i])
	return val

space = 'rspace' #'zspace'
istep = 266

# Path to info
halodir = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/'

# Get the conversion between the name of the time step and redshift                   
step = np.genfromtxt(halodir+'step_redshift.txt',usecols=0,dtype=int)
redshift = np.genfromtxt(halodir+'step_redshift.txt',usecols=1)
zz = redshift[np.where(step == istep)][0]

# Read the HMF mass bins
mfile = halodir + 'hmfs/hmf_'+str(istep)+'.txt'
mf, lmh = np.loadtxt(mfile, usecols= (2, 3), unpack=True)

# Read the bias
mb, bias = np.loadtxt(halodir+'bias_rl20.0_rh80.0.txt', unpack = True)
ibmax = 76
print('Excluding the last {} points'.format(len(mb)-len(mb[:ibmax])))

# Best polynomial fits
p2 = np.poly1d(np.polyfit(mb[:ibmax], bias[:ibmax], 2))
p3 = np.poly1d(np.polyfit(mb[:ibmax], bias[:ibmax], 3))
p4 = np.poly1d(np.polyfit(mb[:ibmax], bias[:ibmax], 4))
p5 = np.poly1d(np.polyfit(mb[:ibmax], bias[:ibmax], 5))

xp = np.linspace(np.min(lmh)-0.5,np.max(lmh)+0.5,100)

# Figure http://matplotlib.org/users/gridspec.html
fig = plt.figure(figsize=(7,9))
gs = gridspec.GridSpec(5,1)
gs.update(wspace=0., hspace=0.)

# Ratio plot
axr = plt.subplot(gs[4,0])
axr.set_xlabel("${\\rm log}M_{\\rm halo} [{\\rm M}_{\odot}/h]$")
#axr.set_ylabel("$b/p_{k}$")
axr.set_ylabel("$b/p_{4}$")
axr.set_autoscale_on(False) ;  axr.minorticks_on()
xmin = 10.6 ; xmax = 14.4
axr.set_xlim(xmin,xmax) ; axr.set_ylim(0.8,1.2)

axr.plot(mb,bias/bias,'k-')
#axr.plot(mb,bias/p2(mb),linestyle='-')
#axr.plot(mb,bias/p3(mb),linestyle='--')
axr.plot(mb,bias/p4(mb),linestyle='-.')
#axr.plot(mb,bias/p5(mb),linestyle=':')

# Bias plot
axb = plt.subplot(gs[2:4,0],sharex=axr)
plt.setp(axb.get_xticklabels(), visible=False)
axb.set_autoscale_on(False) ;  axb.minorticks_on()
axb.set_ylim(0.7,11.)
axb.set_yscale('log')
axb.set_ylabel("$b$")

axb.plot(mb,bias,'ko',label='Haloes')
#axb.plot(xp,p2(xp),label='$p_2$',linestyle='-')
#axb.plot(xp,p3(xp),label='$p_3$',linestyle='--')
axb.plot(xp,p4(xp),label='$p_4$',linestyle='-.')
#axb.plot(xp,p5(xp),label='$p_5$',linestyle=':')

leg = axb.legend(loc=2)
leg.draw_frame(False)

# HMF
axh = plt.subplot(gs[0:2,0],sharex=axb)
plt.setp(axh.get_xticklabels(), visible=False)
axh.set_autoscale_on(False) ;  axh.minorticks_on()
ymin = 0.0000002 ; ymax = 1.
axh.set_ylim(ymin, ymax)
axh.set_yscale('log')
axh.set_ylabel("${\\rm d}n(h^3{\\rm Mpc}^{-3})/{\\rm dlog}M}$")
axh.annotate('z='+str(zz),xy=(13.5,0.1))

ind = np.where(mf >0.)
axh.plot(lmh[ind],mf[ind],'k-')

axr.tick_params(which='both',bottom=True, top=True, left=True, right=True,direction='in')
axb.tick_params(which='both',bottom=True, top=True, left=True, right=True,direction='in')
axh.tick_params(which='both',bottom=True, top=True, left=True, right=True,direction='in')

# Save figure
plotfile = halodir + 'plots/hmf_bias_'+str(istep)+'.png'
fig.savefig(plotfile,dpi=300)
print('Output: ',plotfile)
