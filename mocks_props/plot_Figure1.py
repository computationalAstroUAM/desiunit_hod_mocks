import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
from math import pi
import warnings
import mpl_style
import matplotlib.gridspec as gridspec
plt.style.use(mpl_style.style1)
edges = np.concatenate( [ np.array( np.arange(10.58,12.4999,0.04)), np.array(np.arange(12.5,13.6,0.05)),  np.array(np.arange(13.6,14.0,0.1)),np.array(np.arange(14.2,14.6,0.2)) ])

dm = edges[1:]-edges[:-1]
mhist = edges[1:]-0.5*dm
print(len(mhist))
mass_bins = len(mhist)
bias = np.zeros((mass_bins))

for i in range(0,mass_bins):
    files = str(i)+'.txt'
    bias[i] = np.genfromtxt('../../DESI_outputs/bias_snap100/bias_'+files)


bias = np.delete(bias,76)
mhist = np.delete(mhist,76)

a = np.array([mhist,bias]).T
np.savetxt('../../DESI_outputs/bias_snap100/all_bias.txt',a)

# Best polynomial fits
maxim = len(edges)-2
p2 = np.poly1d(np.polyfit(mhist[:maxim], bias[:maxim], 2))
p3 = np.poly1d(np.polyfit(mhist[:maxim], bias[:maxim], 3))
p4 = np.poly1d(np.polyfit(mhist[:maxim], bias[:maxim], 4))
p5 = np.poly1d(np.polyfit(mhist[:maxim], bias[:maxim], 5))
p6 = np.poly1d(np.polyfit(mhist[:maxim], bias[:maxim], 6))

mb, bias2 = np.loadtxt('bias_rl20.0_rh80.0.txt', unpack = True)

#print('My mass: ',mhist)
#print('Avila mass: ',mb)
#print('My bias: ',bias)
#print('Avila bias: ', bias2)

data = np.loadtxt('../../DESI_outputs/output/HMF_and_accumulated_and_log_mainhalos_snap100.txt')

m = data[:,0]
hmf_log = data[:,3]

fig, axs = plt.subplots(3, 1,sharex=True,gridspec_kw={'height_ratios': [2,4,1]})
fig.subplots_adjust(hspace=0)
axs[0].semilogy(m, hmf_log,'b-',label='UNITSIM',basey=10)
axs[0].set_xlabel(r'$log_{10}\left(M_{h}\right) h^{-1}M_{\odot}$')
axs[0].set_ylabel(r'$\frac{dn}{dlog_{10}M_h} \left(M_{h}\right)$')
axs[0].text(11,10**(-5),r'$z = 0.8594$')
axs[0].legend()

axs[1].semilogy(mhist,bias,'ko',label='Haloes')
axs[1].semilogy(mb, bias2, color='g',label='Avila (2020)')
#axs[1].semilogy(mhist,p2(mhist),color='violet',linestyle='solid',label=r'$p_2$')
#axs[1].semilogy(mhist,p3(mhist),color='g',linestyle='dotted',label=r'$p_3$')
#axs[1].semilogy(mhist,p4(mhist),color='b',linestyle='dashed',label=r'$p_4$')
#axs[1].semilogy(mhist,p5(mhist),color='r',linestyle='dashdot',label=r'$p_5$')
axs[1].semilogy(mhist,p6(mhist),color='y',linestyle='dashed',label=r'$p_6$')
axs[1].set_xlabel(r'$log_{10}M_h \left[h^{-1}M_{\odot}\right]$')
axs[1].set_xlim(10.5,14.5)
axs[1].set_ylabel(r'$b(M_h)$')
axs[1].legend()

axs[2].plot(mhist,bias/bias,'k-')
axs[2].plot(mhist,np.interp(mhist,mb, bias2)/bias, color='green')
#axs[2].plot(mhist,bias/p2(mhist),color='violet',linestyle='solid')
#axs[2].plot(mhist,bias/p3(mhist),color='g',linestyle='dotted')
#axs[2].plot(mhist,bias/p4(mhist),color='b',linestyle='dashed')
#axs[2].plot(mhist,bias/p5(mhist),color='r',linestyle='dashdot')
axs[2].plot(mhist,bias/p6(mhist),color='y',linestyle='dashed')
axs[2].set_xlabel(r'$log_{10}M_h \left[h^{-1}M_{\odot}\right]$')
axs[2].set_xlim(10.5,14.5)
axs[2].set_ylabel(r'$\frac{b}{p_k}$')
plt.savefig('../../DESI_outputs/plots/Figure1_Avila2020_snap100.png',dpi=300)
plt.close()


#fig, axs = plt.subplots(2, 1,sharex=True)
#fig.subplots_adjust(hspace=0)
#axs[0].loglog(10**b, a,'b',label='Mass function',basex=10)
#axs[0].loglog(m2,dndm12,'g--',label='Behroozi reference',basex=10)
#axs[0].loglog(m2,dndm22,'r-.',label='Angulo',basex=10)
#axs[0].loglog(m2,dndm32,'k:',label='Angulo subhaloes',basex=10)
#axs[0].set_ylabel(r'$\frac{dn}{dm}$ $\left[h^{4}Mpc^{-3}M_\odot^{-1}\right]$')
#axs[0].legend()

#axs[1].semilogx(m2, (np.interp(m2,10**b,a)/dndm12),'b',label='Mass function',basex=10)
#axs[1].semilogx(m2,dndm12/dndm12,'g--',label='Behroozi reference',basex=10)
#axs[1].semilogx(m2,dndm22/dndm12,'r-.',label='Angulo',basex=10)
#axs[1].semilogx(m2,dndm32/dndm12,'k:',label='Angulo subhaloes',basex=10)
#axs[1].set_xlabel(r'$M_h$ $\left[h^{-1} M_{\odot}\right]$')
#axs[1].set_ylabel(r'$\frac{\frac{dn}{dm}}{\left(\frac{dn}{dm}\right)_{Behr}}$')
#plt.savefig('plots/Funciondemasa_box1000_comparacion2.png',dpi=300)
#close()
















