import numpy as np
import camb as cb
from camb import get_matter_power_interpolator
from nbodykit.source.catalog import CSVCatalog
from nbodykit.lab import FFTPower #Fast Fourier Transform, tecnica muy estandar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
from math import pi
from scipy.special import spherical_jn #importamos la funcion esferica de bessel que necesitaremos despues
import scipy.integrate as spi
import scipy.integrate as integrate


data1 = np.loadtxt('data_multipoles/n_25e-4_mock1.txt')

k = data1[:,0]
P_0_data1 = data1[:,1]
P_2_data1 = data1[:,2]
P_4_data1 = data1[:,3]

data2 = np.loadtxt('data_multipoles/n_25e-4_mock2.txt')

P_0_data2 = data2[:,1]
P_2_data2 = data2[:,2]
P_4_data2 = data2[:,3]

data3 = np.loadtxt('data_multipoles/n_25e-4_mock3.txt')

P_0_data3 = data3[:,1]
P_2_data3 = data3[:,2]
P_4_data3 = data3[:,3]

data4 = np.loadtxt('data_multipoles/n_25e-4_mock4.txt')

P_0_data4 = data4[:,1]
P_2_data4 = data4[:,2]
P_4_data4 = data4[:,3]

data5 = np.loadtxt('data_multipoles/n_25e-4_mock9.txt')

P_0_data5 = data5[:,1]
P_2_data5 = data5[:,2]
P_4_data5 = data5[:,3]

data6 = np.loadtxt('data_multipoles/n_20e-4_mock1.txt')

P_0_data6 = data6[:,1]
P_2_data6 = data6[:,2]
P_4_data6 = data6[:,3]

data7 = np.loadtxt('data_multipoles/n_20e-4_mock2.txt')

P_0_data7 = data7[:,1]
P_2_data7 = data7[:,2]
P_4_data7 = data7[:,3]

data8 = np.loadtxt('data_multipoles/n_20e-4_mock3.txt')

P_0_data8 = data8[:,1]
P_2_data8 = data8[:,2]
P_4_data8 = data8[:,3]

data9 = np.loadtxt('data_multipoles/n_20e-4_mock4.txt')

P_0_data9 = data9[:,1]
P_2_data9 = data9[:,2]
P_4_data9 = data9[:,3]

data10 = np.loadtxt('data_multipoles/n_20e-4_mock9.txt')

P_0_data10 = data10[:,1]
P_2_data10 = data10[:,2]
P_4_data10 = data10[:,3]

data11 = np.loadtxt('data_multipoles/n_5e-4_mock1.txt')

P_0_data11 = data11[:,1]
P_2_data11 = data11[:,2]
P_4_data11 = data11[:,3]

data12 = np.loadtxt('data_multipoles/n_5e-4_mock2.txt')

P_0_data12 = data12[:,1]
P_2_data12 = data12[:,2]
P_4_data12 = data12[:,3]

data13 = np.loadtxt('data_multipoles/n_5e-4_mock3.txt')

P_0_data13 = data13[:,1]
P_2_data13 = data13[:,2]
P_4_data13 = data13[:,3]

data14 = np.loadtxt('data_multipoles/n_5e-4_mock4.txt')

P_0_data14 = data14[:,1]
P_2_data14 = data14[:,2]
P_4_data14 = data14[:,3]

data15 = np.loadtxt('data_multipoles/n_5e-4_mock9.txt')

P_0_data15 = data15[:,1]
P_2_data15 = data15[:,2]
P_4_data15 = data15[:,3]

sel = np.where((k < 0.55))

k = k[sel]
P_0_data1 = P_0_data1[sel]
P_2_data1 = P_2_data1[sel]
P_4_data1 = P_4_data1[sel]
P_0_data2 = P_0_data2[sel]
P_2_data2 = P_2_data2[sel]
P_4_data2 = P_4_data2[sel]
P_0_data3 = P_0_data3[sel]
P_2_data3 = P_2_data3[sel]
P_4_data3 = P_4_data3[sel]
P_0_data4 = P_0_data4[sel]
P_2_data4 = P_2_data4[sel]
P_4_data4 = P_4_data4[sel]
P_0_data5 = P_0_data5[sel]
P_2_data5 = P_2_data5[sel]
P_4_data5 = P_4_data5[sel]
P_0_data6 = P_0_data6[sel]
P_2_data6 = P_2_data6[sel]
P_4_data6 = P_4_data6[sel]
P_0_data7 = P_0_data7[sel]
P_2_data7 = P_2_data7[sel]
P_4_data7 = P_4_data7[sel]
P_0_data8 = P_0_data8[sel]
P_2_data8 = P_2_data8[sel]
P_4_data8 = P_4_data8[sel]
P_0_data9 = P_0_data9[sel]
P_2_data9 = P_2_data9[sel]
P_4_data9 = P_4_data9[sel]
P_0_data10 = P_0_data10[sel]
P_2_data10 = P_2_data10[sel]
P_4_data10 = P_4_data10[sel]
P_0_data11 = P_0_data11[sel]
P_2_data11 = P_2_data11[sel]
P_4_data11 = P_4_data11[sel]
P_0_data12 = P_0_data12[sel]
P_2_data12 = P_2_data12[sel]
P_4_data12 = P_4_data12[sel]
P_0_data13 = P_0_data13[sel]
P_2_data13 = P_2_data13[sel]
P_4_data13 = P_4_data13[sel]
P_0_data14 = P_0_data14[sel]
P_2_data14 = P_2_data14[sel]
P_4_data14 = P_4_data14[sel]
P_0_data15 = P_0_data15[sel]
P_2_data15 = P_2_data15[sel]
P_4_data15 = P_4_data15[sel]





plt.plot(k,P_0_data1*k,color='b',label='Mock 1')
plt.plot(k,P_0_data2*k,color='y',label='Mock 2')
plt.plot(k,P_0_data3*k,color='r',label='Mock 3')
plt.plot(k,P_0_data4*k,color='g',label='Mock 4')
plt.plot(k,P_0_data5*k,color='k',label='Mock 9')
plt.plot(k,P_2_data1*k,color='b',linestyle ='dashed')
plt.plot(k,P_2_data2*k,color='y',linestyle ='dashed')
plt.plot(k,P_2_data3*k,color='r',linestyle ='dashed')
plt.plot(k,P_2_data4*k,color='g',linestyle ='dashed')
plt.plot(k,P_2_data5*k,color='k',linestyle ='dashed')
plt.plot(k,P_4_data1*k,color='b',linestyle ='dashdot')
plt.plot(k,P_4_data2*k,color='y',linestyle ='dashdot')
plt.plot(k,P_4_data3*k,color='r',linestyle ='dashdot')
plt.plot(k,P_4_data4*k,color='g',linestyle ='dashdot')
plt.plot(k,P_4_data5*k,color='k',linestyle ='dashdot')
plt.xlabel(r'$k$ $\left[h/Mpc\right]$')
plt.ylabel(r'$k\cdot P_l(k)$')
plt.legend()
plt.title('Multipoles n=25e-4')
plt.savefig('plots/P_l_n_25e-4.png',dpi=300)
plt.close()

plt.plot(k,P_0_data6*k,color='b',label='Mock 1')
plt.plot(k,P_0_data7*k,color='y',label='Mock 2')
plt.plot(k,P_0_data8*k,color='r',label='Mock 3')
plt.plot(k,P_0_data9*k,color='g',label='Mock 4')
plt.plot(k,P_0_data10*k,color='k',label='Mock 9')
plt.plot(k,P_2_data6*k,color='b',linestyle ='dashed')
plt.plot(k,P_2_data7*k,color='y',linestyle ='dashed')
plt.plot(k,P_2_data8*k,color='r',linestyle ='dashed')
plt.plot(k,P_2_data9*k,color='g',linestyle ='dashed')
plt.plot(k,P_2_data10*k,color='k',linestyle ='dotted')
plt.plot(k,P_4_data6*k,color='b',linestyle ='dashdot')
plt.plot(k,P_4_data7*k,color='y',linestyle ='dashdot')
plt.plot(k,P_4_data8*k,color='r',linestyle ='dashdot')
plt.plot(k,P_4_data9*k,color='g',linestyle ='dashdot')
plt.plot(k,P_4_data10*k,color='k',linestyle ='dashdot')
plt.xlabel(r'$k$ $\left[h/Mpc\right]$')
plt.ylabel(r'$k\cdot P_l(k)$')
plt.legend()
plt.title('Multipoles n=20e-4')
plt.savefig('plots/P_l_n_20e-4.png',dpi=300)
plt.close()

plt.plot(k,P_0_data11*k,color='b',label='Mock 1')
plt.plot(k,P_0_data12*k,color='y',label='Mock 2')
plt.plot(k,P_0_data13*k,color='r',label='Mock 3')
plt.plot(k,P_0_data14*k,color='g',label='Mock 4')
plt.plot(k,P_0_data15*k,color='k',label='Mock 9')
plt.plot(k,P_2_data11*k,color='b',linestyle ='dashed')
plt.plot(k,P_2_data12*k,color='y',linestyle ='dashed')
plt.plot(k,P_2_data13*k,color='r',linestyle ='dashed')
plt.plot(k,P_2_data14*k,color='g',linestyle ='dashed')
plt.plot(k,P_2_data15*k,color='k',linestyle ='dashed')
plt.plot(k,P_4_data11*k,color='b',linestyle ='dashdot')
plt.plot(k,P_4_data12*k,color='y',linestyle ='dashdot')
plt.plot(k,P_4_data13*k,color='r',linestyle ='dashdot')
plt.plot(k,P_4_data14*k,color='g',linestyle ='dashdot')
plt.plot(k,P_4_data15*k,color='k',linestyle ='dashdot')
plt.xlabel(r'$k$ $\left[h/Mpc\right]$')
plt.ylabel(r'$k\cdot P_l(k)$')
plt.legend()
plt.title('Multipoles n=5e-4')
plt.savefig('plots/P_l_n_5e-4.png',dpi=300)
plt.close()









#Comprobacion fsat=0 en redshift space

data_fsat0 = np.loadtxt('data_multipoles/n_25e-4_mock1_fsat0.txt')

P_0_fsat0 = data_fsat0[:,1]
P_2_fsat0 = data_fsat0[:,2]
P_4_fsat0 = data_fsat0[:,3]

data_fsat0_realspace= np.loadtxt('data_multipoles/n_25e-4_mock1_fsat0_realspace.txt')

P_0_fsat0_realspace = data_fsat0_realspace[:,1]
P_2_fsat0_realspace = data_fsat0_realspace[:,2]
P_4_fsat0_realspace = data_fsat0_realspace[:,3]



DATA_fsat0 =  np.loadtxt('Kaiser_multipoles/n_25e-4_mock1_fsat0_KAISER.txt')

P_0_FSAT0 = DATA_fsat0[:,1]
P_2_FSAT0 = DATA_fsat0[:,2]
P_4_FSAT0 = DATA_fsat0[:,3]

DATA_fsat0_realspace =  np.loadtxt('Kaiser_multipoles/n_25e-4_mock1_fsat0_KAISER_realspace.txt')

P_0_FSAT0_realspace = DATA_fsat0_realspace[:,1]
P_2_FSAT0_realspace = DATA_fsat0_realspace[:,2]
P_4_FSAT0_realspace = DATA_fsat0_realspace[:,3]




P_0_fsat0 = P_0_fsat0[sel]
P_2_fsat0 = P_2_fsat0[sel]
P_4_fsat0 = P_4_fsat0[sel]
P_0_FSAT0 = P_0_FSAT0[sel]
P_2_FSAT0 = P_2_FSAT0[sel]
P_4_FSAT0 = P_4_FSAT0[sel]
P_0_fsat0_realspace = P_0_fsat0_realspace[sel]
P_2_fsat0_realspace = P_2_fsat0_realspace[sel]
P_4_fsat0_realspace = P_4_fsat0_realspace[sel]
P_0_FSAT0_realspace = P_0_FSAT0_realspace[sel]
P_2_FSAT0_realspace = P_2_FSAT0_realspace[sel]
P_4_FSAT0_realspace = P_4_FSAT0_realspace[sel]




plt.semilogx(k,P_0_fsat0/P_0_FSAT0,color='g',linestyle ='dashdot',label=r'$P_0$ $s$')
plt.semilogx(k,P_2_fsat0/P_2_FSAT0,color='r',linestyle ='dashdot',label=r'$P_2$ $s$')
plt.semilogx(k,P_0_fsat0_realspace/P_0_FSAT0_realspace,color='b',linestyle ='solid',label=r'$P_0$ $r$')
#plt.semilogx(k,P_2_fsat0_realspace/P_2_FSAT0_realspace,color='violet',linestyle ='solid',label=r'$P_2$ $r$')
#plt.semilogx(k,P_4_fsat0/P_4_FSAT0,color='b',linestyle ='dashdot',label=r'$P_4$')
plt.axhline(y=1,color='k')
plt.xlabel(r'$k$ $\left[h/Mpc\right]$')
plt.ylabel(r'$P_l (k)/f(b,f)\cdot P_m(k)$')
plt.legend()
plt.title('fsat=0 n=25e-4 multipoles Kaiser')
plt.savefig('plots/P_l_n_25e-4_kaiser_fsat0.png',dpi=300)
plt.close()

realspace_biasrsd = np.loadtxt('Kaiser_multipoles/n_25e-4_mock1_fsat0_KAISER_realspace_biasrsd.txt')

P_0_realspace_biasrsd = realspace_biasrsd[:,1]
P_0_realspace_biasrsd = P_0_realspace_biasrsd[sel]

plt.semilogx(k,P_0_fsat0/P_0_FSAT0,color='g',linestyle ='dashdot',label=r'$P_0$ $s$')
plt.semilogx(k,P_2_fsat0/P_2_FSAT0,color='r',linestyle ='dashdot',label=r'$P_2$ $s$')
plt.semilogx(k,P_0_fsat0_realspace/P_0_realspace_biasrsd,color='b',linestyle ='solid',label=r'$P_0$ $r$')
#plt.semilogx(k,P_2_fsat0_realspace/P_2_FSAT0_realspace,color='violet',linestyle ='solid',label=r'$P_2$ $r$')
#plt.semilogx(k,P_4_fsat0/P_4_FSAT0,color='b',linestyle ='dashdot',label=r'$P_4$')
plt.axhline(y=1,color='k')
plt.xlabel(r'$k$ $\left[h/Mpc\right]$')
plt.ylabel(r'$P_l (k)/f(b,f)\cdot P_m(k)$')
plt.legend()
plt.title('fsat=0 n=25e-4 multipoles Kaiser')
plt.savefig('plots/P_l_n_25e-4_kaiser_fsat0_biasrsd.png',dpi=300)
plt.close()

redshiftspace_bias = np.loadtxt('Kaiser_multipoles/n_25e-4_mock1_fsat0_KAISER_biasrealspace.txt')

P_0_redshiftspace_bias = redshiftspace_bias[:,1]
P_2_redshiftspace_bias = redshiftspace_bias[:,2]

P_0_redshiftspace_bias = P_0_redshiftspace_bias[sel]
P_2_redshiftspace_bias = P_2_redshiftspace_bias[sel]

plt.semilogx(k,P_0_fsat0/P_0_redshiftspace_bias,color='g',linestyle ='dashdot',label=r'$P_0$ $s$')
plt.semilogx(k,P_2_fsat0/P_2_redshiftspace_bias,color='r',linestyle ='dashdot',label=r'$P_2$ $s$')
plt.semilogx(k,P_0_fsat0_realspace/P_0_FSAT0_realspace,color='b',linestyle ='solid',label=r'$P_0$ $r$')
#plt.semilogx(k,P_2_fsat0_realspace/P_2_FSAT0_realspace,color='violet',linestyle ='solid',label=r'$P_2$ $r$')
#plt.semilogx(k,P_4_fsat0/P_4_FSAT0,color='b',linestyle ='dashdot',label=r'$P_4$')
plt.axhline(y=1,color='k')
plt.xlabel(r'$k$ $\left[h/Mpc\right]$')
plt.ylabel(r'$P_l (k)/f(b,f)\cdot P_m(k)$')
plt.legend()
plt.title('fsat=0 n=25e-4 multipoles Kaiser')
plt.savefig('plots/P_l_n_25e-4_kaiser_fsat0_bias.png',dpi=300)
plt.close()









DATA1 = np.loadtxt('Kaiser_multipoles/n_25e-4_mock1_KAISER.txt')

P_0_DATA1 = DATA1[:,1]
P_2_DATA1 = DATA1[:,2]

DATA2 = np.loadtxt('Kaiser_multipoles/n_25e-4_mock2_KAISER.txt')

P_0_DATA2 = DATA2[:,1]
P_2_DATA2 = DATA2[:,2]

DATA3 = np.loadtxt('Kaiser_multipoles/n_25e-4_mock3_KAISER.txt')

P_0_DATA3 = DATA3[:,1]
P_2_DATA3 = DATA3[:,2]

DATA4 = np.loadtxt('Kaiser_multipoles/n_25e-4_mock4_KAISER.txt')

P_0_DATA4 = DATA4[:,1]
P_2_DATA4 = DATA4[:,2]

DATA5 = np.loadtxt('Kaiser_multipoles/n_25e-4_mock9_KAISER.txt')

P_0_DATA5 = DATA5[:,1]
P_2_DATA5 = DATA5[:,2]


P_0_DATA1 = P_0_DATA1[sel]
P_0_DATA2 = P_0_DATA2[sel]
P_0_DATA3 = P_0_DATA3[sel]
P_0_DATA4 = P_0_DATA4[sel]
P_0_DATA5 = P_0_DATA5[sel]
P_2_DATA1 = P_2_DATA1[sel]
P_2_DATA2 = P_2_DATA2[sel]
P_2_DATA3 = P_2_DATA3[sel]
P_2_DATA4 = P_2_DATA4[sel]
P_2_DATA5 = P_2_DATA5[sel]


plt.semilogx(k,P_0_data1/P_0_DATA1,color='b',linestyle ='solid',label='Mock 1')
plt.semilogx(k,P_0_data2/P_0_DATA2,color='y',linestyle ='solid',label='Mock 2')
plt.semilogx(k,P_0_data3/P_0_DATA3,color='r',linestyle ='solid',label='Mock 3')
plt.semilogx(k,P_0_data4/P_0_DATA4,color='g',linestyle ='solid',label='Mock 4')
plt.semilogx(k,P_0_data5/P_0_DATA5,color='k',linestyle ='solid',label='Mock 9')
plt.xlabel(r'$k$ $\left[h/Mpc\right]$')
plt.ylabel(r'$P_0/\left(b^2+\frac{2}{3} f\cdot b + \frac{1}{5} f^2\right)\cdot P_m(k)$')
plt.legend()
plt.title('Kaiser effect monopole n=25e-4')
plt.savefig('plots/Kaiser_monopole_n_25e-4.png',dpi=300)
plt.close()

plt.loglog(k,P_0_data1,color='b',linestyle ='solid',label='Mock 1')
plt.loglog(k,P_0_data2,color='y',linestyle ='solid',label='Mock 2')
plt.loglog(k,P_0_data3,color='r',linestyle ='solid',label='Mock 3')
plt.loglog(k,P_0_data4,color='g',linestyle ='solid',label='Mock 4')
plt.loglog(k,P_0_data5,color='k',linestyle ='solid',label='Mock 9')
plt.loglog(k,P_0_DATA1,color='b',linestyle ='dashed',label='Th 1')
plt.loglog(k,P_0_DATA2,color='y',linestyle ='dashed',label='Th 2')
plt.loglog(k,P_0_DATA3,color='r',linestyle ='dashed',label='Th 3')
plt.loglog(k,P_0_DATA4,color='g',linestyle ='dashed',label='Th 4')
plt.loglog(k,P_0_DATA5,color='k',linestyle ='dashed',label='Th 9')
plt.xlabel(r'$k$ $\left[h/Mpc\right]$')
plt.ylabel(r'$P_0/\left(b^2+\frac{2}{3} f\cdot b + \frac{1}{5} f^2\right)\cdot P_m(k)$')
plt.legend()
plt.title('Kaiser effect monopole n=25e-4')
plt.savefig('plots/Kaiser_monopole_loglog_n_25e-4.png',dpi=300)
plt.close()






