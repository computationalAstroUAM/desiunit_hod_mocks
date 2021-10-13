import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
from math import pi
import mpl_style
plt.style.use(mpl_style.style1)

data = np.loadtxt('../../DESI_outputs/HODs/HOD_gaussPL_PL_corr_25e-4.params')
#data = np.loadtxt('fit_mocks/HODs/HOD_gaussPL_PL_corr_20e-4.params')
#data = np.loadtxt('fit_mocks/HODs/HOD_gaussPL_PL_corr_5e-4.params')

f_sat = data[:,0]
mu = data[:,1]
A_c = data[:,2]
A_s = data[:,3]


fig, axs = plt.subplots(3, 1,sharex=True)
fig.subplots_adjust(hspace=0)
axs[0].plot(f_sat, mu,'b-')
axs[0].set_xlabel(r'$f_{sat}$')
axs[0].set_ylabel(r'$\mu$')
axs[0].text(0.6,11.5,r'$z = 0.9873$')

axs[1].plot(f_sat,A_c,'r-')
axs[1].set_xlabel(r'$f_{sat}$')
axs[1].set_ylabel(r'$A_c$')

axs[2].plot(f_sat,A_s,'g-')
axs[2].set_xlabel(r'$f_{sat}$')
axs[2].set_ylabel(r'$A_s$')
#plt.savefig('plots/Ac_As_mu.png',dpi=300)
plt.close()

sel1= np.where((f_sat == 0.21))
sel2= np.where((f_sat == 0.48))
sel3= np.where((f_sat == 0.51))
sel4= np.where((f_sat == 0.56))

mu1 = float(mu[sel1])
A_c1 = float(A_c[sel1])
A_s1 = float(A_s[sel1])
mu2 = float(mu[sel2])
A_c2 = float(A_c[sel2])
A_s2 = float(A_s[sel2])
mu3 = float(mu[sel3])
A_c3 = float(A_c[sel3])
A_s3 = float(A_s[sel3])
mu4 = float(mu[sel4])
A_c4 = float(A_c[sel4])
A_s4 = float(A_s[sel4])

edges = np.concatenate( [ np.array( np.arange(10.58,12.4999,0.04)), np.array(np.arange(12.5,13.6,0.05)),  np.array(np.arange(13.6,14.0,0.1)),np.array(np.arange(14.2,14.6,0.2))])
dm = edges[1:]-edges[:-1]
mhist = edges[1:]-0.5*dm

M = 10**(mhist[:76])
alpha = 0.9  #HOD3
sigma = 0.12  #HOD3
gamma = -1.4   #HOD3

M_0_1 = float(10**(mu1-0.05)) #HOD3
M_1_1 = float(10**(mu1+0.35)) #HOD3
M_0_2 = float(10**(mu2-0.05)) #HOD3
M_1_2 = float(10**(mu2+0.35)) #HOD3
M_0_3 = float(10**(mu3-0.05)) #HOD3
M_1_3 = float(10**(mu3+0.35)) #HOD3
M_0_4 = float(10**(mu4-0.05)) #HOD3
M_1_4 = float(10**(mu4+0.35)) #HOD3

print('alpha=',alpha)
print('sigma=',sigma)
print('gamma=',gamma)


N_sat_M_1 = np.zeros((len(M)))
N_sat_M_2 = np.zeros((len(M)))
N_sat_M_3 = np.zeros((len(M)))
N_sat_M_4 = np.zeros((len(M)))

for i, Mi in enumerate(M):
    if np.log10(Mi) > np.log10(M_0_1):
        N_sat_M_1[i] = A_s1*((Mi-M_0_1)/M_1_1)**alpha
    else:
        N_sat_M_1[i] = 0

for i, Mi in enumerate(M):
    if np.log10(Mi) > np.log10(M_0_2):
        N_sat_M_2[i] = A_s2*((Mi-M_0_2)/M_1_2)**alpha
    else:
        N_sat_M_2[i] = 0

for i, Mi in enumerate(M):
    if np.log10(Mi) > np.log10(M_0_3):
        N_sat_M_3[i] = A_s3*((Mi-M_0_3)/M_1_3)**alpha
    else:
        N_sat_M_3[i] = 0

for i, Mi in enumerate(M):
    if np.log10(Mi) > np.log10(M_0_4):
        N_sat_M_4[i] = A_s4*((Mi-M_0_4)/M_1_4)**alpha
    else:
        N_sat_M_4[i] = 0



N_cen_M_1 = np.zeros((len(M)))
N_cen_M_2 = np.zeros((len(M)))
N_cen_M_3 = np.zeros((len(M)))
N_cen_M_4 = np.zeros((len(M)))

const1=A_c1/(np.sqrt(2*np.pi)*sigma)
const2=A_c2/(np.sqrt(2*np.pi)*sigma)
const3=A_c3/(np.sqrt(2*np.pi)*sigma)
const4=A_c4/(np.sqrt(2*np.pi)*sigma)


i=0
for i,Mi in enumerate(M):
    if np.log10(Mi) < mu1:
        N_cen_M_1[i] = const1*np.exp(-((np.log10(Mi)-mu1)**2)/(2*sigma**2))
    else:
        N_cen_M_1[i] = const1*(Mi/(10**mu1))**gamma

for i,Mi in enumerate(M):
    if np.log10(Mi) < mu2:
        N_cen_M_2[i] = const2*np.exp(-((np.log10(Mi)-mu2)**2)/(2*sigma**2))
    else:
        N_cen_M_2[i] = const2*(Mi/(10**mu2))**gamma

for i,Mi in enumerate(M):
    if np.log10(Mi) < mu3:
        N_cen_M_3[i] = const3*np.exp(-((np.log10(Mi)-mu3)**2)/(2*sigma**2))
    else:
        N_cen_M_3[i] = const3*(Mi/(10**mu3))**gamma

for i,Mi in enumerate(M):
    if np.log10(Mi) < mu4:
        N_cen_M_4[i] = const4*np.exp(-((np.log10(Mi)-mu4)**2)/(2*sigma**2))
    else:
        N_cen_M_4[i] = const4*(Mi/(10**mu4))**gamma




plt.loglog(M,N_sat_M_1,'b-.',label=r'$f_{sat}=0.21$ $HOD3$')
plt.loglog(M,N_cen_M_1,'b.')
plt.loglog(M,N_sat_M_1+N_cen_M_1,'b-')
plt.loglog(M,N_sat_M_2,'y-.',label=r'$f_{sat}=0.48$ $HOD3$')
plt.loglog(M,N_cen_M_2,'y.')
plt.loglog(M,N_sat_M_2+N_cen_M_2,'y-')
plt.loglog(M,N_sat_M_3,'g-.',label=r'$f_{sat}=0.51$ $HOD3$')
plt.loglog(M,N_cen_M_3,'g.')
plt.loglog(M,N_sat_M_3+N_cen_M_3,'g-')
plt.loglog(M,N_sat_M_4,'r-.',label=r'$f_{sat}=0.56$ $HOD3$')
plt.loglog(M,N_cen_M_4,'r.')
plt.loglog(M,N_sat_M_4+N_cen_M_4,'r-')
plt.xlabel(r'$log_{10}(M)$ $\left[M_\odot/h\right]$')
plt.ylabel(r'$<N>$')
plt.ylim(10**(-6),100)
plt.text(10**13.4,10**(-0.4),r'$\alpha = 0.9$')
plt.text(10**13.4,10**(-0.8),r'$\sigma = 0.12$')
plt.text(10**13.4,10**(-1.2),r'$\gamma = -1.4$')
plt.text(10**13.0,10**(-1.6),r'$b_{gal}=1.4$')
plt.text(10**13.0,10**(-2.0),r'$n_{gal}=25\cdot 10^{-4}$')
plt.legend()
plt.savefig('../../DESI_outputs/plots/N_sat_cen.png',dpi=300)




