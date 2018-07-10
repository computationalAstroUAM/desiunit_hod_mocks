import numpy as np
import scipy.interpolate as interp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
matplotlib.rcParams.update({'font.size': 12})

fig3 = plt.figure(figsize = (8,6))

DM = "/users/savila/CorrelationTools/xir/OuterRim.b1.00.z0.865.2PCF"
DMarr = np.loadtxt(DM)

fig = plt.figure(figsize = (8,6))
ax = fig.add_subplot(111)
#FILENAME= "../output/HOD_V3/galaxies_V1.0_M1_12.00_Mmin11.00_M0_12.00_alpha_1.00_sig_0.25.2PCF"
FILENAME= "../../../data/Favole15/Combined_W134_xi.CORR.dat"
array= np.loadtxt(FILENAME)
FILENAME= "../eBOSS_LSS_ELGs/xi0gebossELG_eboss211_mz0.7xz1.1fkp8st0.dat"
array= np.loadtxt(FILENAME)
sel= np.array(np.where(array[:,1]>=0))
array = array[sel[0,:]]
#err2PCF = (array[:,1]+1.)/(np.sqrt(array[:,2])*array[:,1])
DMinterp = np.interp(array[:,0], DMarr[:,0], DMarr[:,1])
b = np.sqrt(array[:,1]/DMinterp)
#errb = np.sqrt(err2PCF*err2PCF/(4.*DMinterp*array[:,1]))
#w = 1.0/(errb[np.where(array[:,0]>5.0)]**2) * np.sum(errb[np.where(array[:,0]>5.0)]**2)
#linfit=np.polyfit(array[np.where(array[:,0]>5.0)[0], 0], b[np.where(array[:,0]>5.0)],deg=0)#, w=w)
		#print('DLM=',DLM[i],'M1=',M1[j],'mean b =',linfit)
#		b_10 = np.interp(10.0,array[:,0], b)
#ax.plot(array[:,0], b , marker="o" , markersize=3.0 ,color=colors[j],linestyle = 'None')
#ax.errorbar(array[:,0],b, yerr=errb)
#		ax.plot(array[:, 0], np.zeros_like(array[:,0])+b_10 ,color=colors[j],linestyle = 'dashed', label='b={}'.format(b_10))
ax.plot(array[:,0], DMinterp)
ax.plot(array[:,0], array)
	#	ax.plot(array[:, 0], np.zeros_like(array[:,0])+mean ,color=colors[j],linestyle = 'dashed', label='b={}'.format(mean))
#	ax.plot(HALOarr[sel2[0,:],0], HALOarr[sel2[0,:],1],"k-", label='Halos')
#	ax.plot(eBOSSarr[sel3[0,:],0], eBOSSarr[sel3[0,:],1],"kx",markersize=5, label='eBOSS')
#ax.set_xscale("log")
#ax.set_yscale("log")
ax.set_ylim(1, 3)
ax.legend(fontsize=7,loc='upper right' ) 
ax.set_xlabel(r"Comoving Distance r[$Mpc h^{-1}$]")
ax.set_ylabel(r"Galaxy Bias $b(r)$")
plt.savefig("HOD_V3/GalBias_Favole.png")


    
