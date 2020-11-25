import numpy as np
import sys
from math import pi

Lbox = 1000 #Mpc/h
H_0 = 100 #Mpc/h
Om = 0.3089

def H(z):
    return H_0*np.sqrt(Om*(1+z)**3+1-Om)

z = 0.9873


n_25_mock1 = np.loadtxt('hod_mocks/output_V1/n_25e-4/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0089_As0.01954_vfact1.00_beta-2.000_K1.00_vt0pm0_mock.dat')

X_25_1 = n_25_mock1[:,0]
Y_25_1 = n_25_mock1[:,1]
Z_25_1 = n_25_mock1[:,2]
VZ_25_1 = n_25_mock1[:,5]
Z_RSD_25_1 = X_25_1+(1+z)*VZ_25_1/H(z)

for i in range(0,len(Z_RSD_25_1)):
    if Z_RSD_25_1[i] < 0:
        Z_RSD_25_1[i] = Z_RSD_25_1[i]+Lbox
    if Z_RSD_25_1[i] > Lbox:
        Z_RSD_25_1[i] = Z_RSD_25_1[i]-Lbox

out_25_1 = np.array([X_25_1,Y_25_1,Z_25_1,Z_RSD_25_1]).T
np.savetxt('NERSC_tables/n_25_mock1.txt',out_25_1)



n_25_mock2 = np.loadtxt('hod_mocks/output_V1/n_25e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.079_Ac0.0129_As0.02471_vfact1.00_beta0.000_K0.25_vt0pm0_mock.dat')

X_25_2 = n_25_mock2[:,0]
Y_25_2 = n_25_mock2[:,1]
Z_25_2 = n_25_mock2[:,2]
VZ_25_2 = n_25_mock2[:,5]
Z_RSD_25_2 = X_25_2+(1+z)*VZ_25_2/H(z)

for i in range(0,len(Z_RSD_25_2)):
    if Z_RSD_25_2[i] < 0:
        Z_RSD_25_2[i] = Z_RSD_25_2[i]+Lbox
    if Z_RSD_25_2[i] > Lbox:
        Z_RSD_25_2[i] = Z_RSD_25_2[i]-Lbox

out_25_2 = np.array([X_25_2,Y_25_2,Z_25_2,Z_RSD_25_2]).T
np.savetxt('NERSC_tables/n_25_mock2.txt',out_25_2)



n_25_mock3 = np.loadtxt('hod_mocks/output_V1/n_25e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.723_Ac0.0860_As0.05916_vfact1.50_beta0.000_K1.00_vt0pm0_mock.dat')

X_25_3 = n_25_mock3[:,0]
Y_25_3 = n_25_mock3[:,1]
Z_25_3 = n_25_mock3[:,2]
VZ_25_3 = n_25_mock3[:,5]
Z_RSD_25_3 = X_25_3+(1+z)*VZ_25_3/H(z)

for i in range(0,len(Z_RSD_25_3)):
    if Z_RSD_25_3[i] < 0:
        Z_RSD_25_3[i] = Z_RSD_25_3[i]+Lbox
    if Z_RSD_25_3[i] > Lbox:
        Z_RSD_25_3[i] = Z_RSD_25_3[i]-Lbox

out_25_3 = np.array([X_25_3,Y_25_3,Z_25_3,Z_RSD_25_3]).T
np.savetxt('NERSC_tables/n_25_mock3.txt',out_25_3)



n_25_mock4 = np.loadtxt('hod_mocks/output_V1/n_25e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.723_Ac0.0860_As0.05916_vfact1.00_beta0.000_K1.00_vt500pm200_mock.dat')

X_25_4 = n_25_mock4[:,0]
Y_25_4 = n_25_mock4[:,1]
Z_25_4 = n_25_mock4[:,2]
VZ_25_4 = n_25_mock4[:,5]
Z_RSD_25_4 = X_25_4+(1+z)*VZ_25_4/H(z)

for i in range(0,len(Z_RSD_25_4)):
    if Z_RSD_25_4[i] < 0:
        Z_RSD_25_4[i] = Z_RSD_25_4[i]+Lbox
    if Z_RSD_25_4[i] > Lbox:
        Z_RSD_25_4[i] = Z_RSD_25_4[i]-Lbox

out_25_4 = np.array([X_25_4,Y_25_4,Z_25_4,Z_RSD_25_4]).T
np.savetxt('NERSC_tables/n_25_mock4.txt',out_25_4)



n_25_mock9 = np.loadtxt('hod_mocks/output_V1/n_25e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.153_Ac0.0160_As0.02829_vfact1.00_beta0.100_K0.15_vt0pm0_mock.dat')

X_25_5 = n_25_mock9[:,0]
Y_25_5 = n_25_mock9[:,1]
Z_25_5 = n_25_mock9[:,2]
VZ_25_5 = n_25_mock9[:,5]
Z_RSD_25_5 = X_25_5+(1+z)*VZ_25_5/H(z)

for i in range(0,len(Z_RSD_25_5)):
    if Z_RSD_25_5[i] < 0:
        Z_RSD_25_5[i] = Z_RSD_25_5[i]+Lbox
    if Z_RSD_25_5[i] > Lbox:
        Z_RSD_25_5[i] = Z_RSD_25_5[i]-Lbox

out_25_5 = np.array([X_25_5,Y_25_5,Z_25_5,Z_RSD_25_5]).T
np.savetxt('NERSC_tables/n_25_mock9.txt',out_25_5)




















n_20_mock1 = np.loadtxt('hod_mocks/output_V1/n_20e-4/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0071_As0.01563_vfact1.00_beta-2.000_K1.00_vt0pm0_mock.dat')

X_20_1 = n_20_mock1[:,0]
Y_20_1 = n_20_mock1[:,1]
Z_20_1 = n_20_mock1[:,2]
VZ_20_1 = n_20_mock1[:,5]
Z_RSD_20_1 = X_20_1+(1+z)*VZ_20_1/H(z)

for i in range(0,len(Z_RSD_20_1)):
    if Z_RSD_20_1[i] < 0:
        Z_RSD_20_1[i] = Z_RSD_20_1[i]+Lbox
    if Z_RSD_20_1[i] > Lbox:
        Z_RSD_20_1[i] = Z_RSD_20_1[i]-Lbox

out_20_1 = np.array([X_20_1,Y_20_1,Z_20_1,Z_RSD_20_1]).T
np.savetxt('NERSC_tables/n_20_mock1.txt',out_20_1)



n_20_mock2 = np.loadtxt('hod_mocks/output_V1/n_20e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.079_Ac0.0103_As0.01977_vfact1.00_beta0.000_K0.25_vt0pm0_mock.dat')

X_20_2 = n_20_mock2[:,0]
Y_20_2 = n_20_mock2[:,1]
Z_20_2 = n_20_mock2[:,2]
VZ_20_2 = n_20_mock2[:,5]
Z_RSD_20_2 = X_20_2+(1+z)*VZ_20_2/H(z)

for i in range(0,len(Z_RSD_20_2)):
    if Z_RSD_20_2[i] < 0:
        Z_RSD_20_2[i] = Z_RSD_20_2[i]+Lbox
    if Z_RSD_20_2[i] > Lbox:
        Z_RSD_20_2[i] = Z_RSD_20_2[i]-Lbox

out_20_2 = np.array([X_20_2,Y_20_2,Z_20_2,Z_RSD_20_2]).T
np.savetxt('NERSC_tables/n_20_mock2.txt',out_20_2)



n_20_mock3 = np.loadtxt('hod_mocks/output_V1/n_20e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.723_Ac0.0688_As0.04733_vfact1.50_beta0.000_K1.00_vt0pm0_mock.dat')

X_20_3 = n_20_mock3[:,0]
Y_20_3 = n_20_mock3[:,1]
Z_20_3 = n_20_mock3[:,2]
VZ_20_3 = n_20_mock3[:,5]
Z_RSD_20_3 = X_20_3+(1+z)*VZ_20_3/H(z)

for i in range(0,len(Z_RSD_20_3)):
    if Z_RSD_20_3[i] < 0:
        Z_RSD_20_3[i] = Z_RSD_20_3[i]+Lbox
    if Z_RSD_20_3[i] > Lbox:
        Z_RSD_20_3[i] = Z_RSD_20_3[i]-Lbox

out_20_3 = np.array([X_20_3,Y_20_3,Z_20_3,Z_RSD_20_3]).T
np.savetxt('NERSC_tables/n_20_mock3.txt',out_20_3)




n_20_mock4 = np.loadtxt('hod_mocks/output_V1/n_20e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.723_Ac0.0688_As0.04733_vfact1.00_beta0.000_K1.00_vt500pm200_mock.dat')

X_20_4 = n_20_mock4[:,0]
Y_20_4 = n_20_mock4[:,1]
Z_20_4 = n_20_mock4[:,2]
VZ_20_4 = n_20_mock4[:,5]
Z_RSD_20_4 = X_20_4+(1+z)*VZ_20_4/H(z)

for i in range(0,len(Z_RSD_20_4)):
    if Z_RSD_20_4[i] < 0:
        Z_RSD_20_4[i] = Z_RSD_20_4[i]+Lbox
    if Z_RSD_20_4[i] > Lbox:
        Z_RSD_20_4[i] = Z_RSD_20_4[i]-Lbox

out_20_4 = np.array([X_20_4,Y_20_4,Z_20_4,Z_RSD_20_4]).T
np.savetxt('NERSC_tables/n_20_mock4.txt',out_20_4)



n_20_mock9 = np.loadtxt('hod_mocks/output_V1/n_20e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.153_Ac0.0128_As0.02263_vfact1.00_beta0.100_K0.15_vt0pm0_mock.dat')

X_20_5 = n_20_mock9[:,0]
Y_20_5 = n_20_mock9[:,1]
Z_20_5 = n_20_mock9[:,2]
VZ_20_5 = n_20_mock9[:,5]
Z_RSD_20_5 = X_20_5+(1+z)*VZ_20_5/H(z)

for i in range(0,len(Z_RSD_20_5)):
    if Z_RSD_20_5[i] < 0:
        Z_RSD_20_5[i] = Z_RSD_20_5[i]+Lbox
    if Z_RSD_20_5[i] > Lbox:
        Z_RSD_20_5[i] = Z_RSD_20_5[i]-Lbox

out_20_5 = np.array([X_20_5,Y_20_5,Z_20_5,Z_RSD_20_5]).T
np.savetxt('NERSC_tables/n_20_mock9.txt',out_20_5)





n_5_mock1 = np.loadtxt('hod_mocks/output_V1/n_5e-4/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0018_As0.00391_vfact1.00_beta-2.000_K1.00_vt0pm0_mock.dat')

X_5_1 = n_5_mock1[:,0]
Y_5_1 = n_5_mock1[:,1]
Z_5_1 = n_5_mock1[:,2]
VZ_5_1 = n_5_mock1[:,5]
Z_RSD_5_1 = X_5_1+(1+z)*VZ_5_1/H(z)

for i in range(0,len(Z_RSD_5_1)):
    if Z_RSD_5_1[i] < 0:
        Z_RSD_5_1[i] = Z_RSD_5_1[i]+Lbox
    if Z_RSD_5_1[i] > Lbox:
        Z_RSD_5_1[i] = Z_RSD_5_1[i]-Lbox

out_5_1 = np.array([X_5_1,Y_5_1,Z_5_1,Z_RSD_5_1]).T
np.savetxt('NERSC_tables/n_5_mock1.txt',out_5_1)



n_5_mock2 = np.loadtxt('hod_mocks/output_V1/n_5e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.079_Ac0.0026_As0.00494_vfact1.00_beta0.000_K0.25_vt0pm0_mock.dat')

X_5_2 = n_5_mock2[:,0]
Y_5_2 = n_5_mock2[:,1]
Z_5_2 = n_5_mock2[:,2]
VZ_5_2 = n_5_mock2[:,5]
Z_RSD_5_2 = X_5_2+(1+z)*VZ_5_2/H(z)

for i in range(0,len(Z_RSD_5_2)):
    if Z_RSD_5_2[i] < 0:
        Z_RSD_5_2[i] = Z_RSD_5_2[i]+Lbox
    if Z_RSD_5_2[i] > Lbox:
        Z_RSD_5_2[i] = Z_RSD_5_2[i]-Lbox

out_5_2 = np.array([X_5_2,Y_5_2,Z_5_2,Z_RSD_5_2]).T
np.savetxt('NERSC_tables/n_5_mock2.txt',out_5_2)



n_5_mock3 = np.loadtxt('hod_mocks/output_V1/n_5e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.723_Ac0.0172_As0.01183_vfact1.50_beta0.000_K1.00_vt0pm0_mock.dat')

X_5_3 = n_5_mock3[:,0]
Y_5_3 = n_5_mock3[:,1]
Z_5_3 = n_5_mock3[:,2]
VZ_5_3 = n_5_mock3[:,5]
Z_RSD_5_3 = X_5_3+(1+z)*VZ_5_3/H(z)

for i in range(0,len(Z_RSD_5_3)):
    if Z_RSD_5_3[i] < 0:
        Z_RSD_5_3[i] = Z_RSD_5_3[i]+Lbox
    if Z_RSD_5_3[i] > Lbox:
        Z_RSD_5_3[i] = Z_RSD_5_3[i]-Lbox

out_5_3 = np.array([X_5_3,Y_5_3,Z_5_3,Z_RSD_5_3]).T
np.savetxt('NERSC_tables/n_5_mock3.txt',out_5_3)



n_5_mock4 = np.loadtxt('hod_mocks/output_V1/n_5e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.723_Ac0.0172_As0.01183_vfact1.00_beta0.000_K1.00_vt500pm200_mock.dat')

X_5_4 = n_5_mock4[:,0]
Y_5_4 = n_5_mock4[:,1]
Z_5_4 = n_5_mock4[:,2]
VZ_5_4 = n_5_mock4[:,5]
Z_RSD_5_4 = X_5_4+(1+z)*VZ_5_4/H(z)

for i in range(0,len(Z_RSD_5_4)):
    if Z_RSD_5_4[i] < 0:
        Z_RSD_5_4[i] = Z_RSD_5_4[i]+Lbox
    if Z_RSD_5_4[i] > Lbox:
        Z_RSD_5_4[i] = Z_RSD_5_4[i]-Lbox

out_5_4 = np.array([X_5_4,Y_5_4,Z_5_4,Z_RSD_5_4]).T
np.savetxt('NERSC_tables/n_5_mock4.txt',out_5_4)



n_5_mock9 = np.loadtxt('hod_mocks/output_V1/n_5e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.153_Ac0.0032_As0.00566_vfact1.00_beta0.100_K0.15_vt0pm0_mock.dat')

X_5_5 = n_5_mock9[:,0]
Y_5_5 = n_5_mock9[:,1]
Z_5_5 = n_5_mock9[:,2]
VZ_5_5 = n_5_mock9[:,5]
Z_RSD_5_5 = X_5_5+(1+z)*VZ_5_5/H(z)

for i in range(0,len(Z_RSD_5_5)):
    if Z_RSD_5_5[i] < 0:
        Z_RSD_5_5[i] = Z_RSD_5_5[i]+Lbox
    if Z_RSD_5_5[i] > Lbox:
        Z_RSD_5_5[i] = Z_RSD_5_5[i]-Lbox

out_5_5 = np.array([X_5_5,Y_5_5,Z_5_5,Z_RSD_5_5]).T
np.savetxt('NERSC_tables/n_5_mock9.txt',out_5_5)


