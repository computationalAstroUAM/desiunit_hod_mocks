import numpy as np
import sys
from math import pi

H_0 = 100 #Mpc/h

#./HOD_NFW_V14_more_BVG 10.944 0.00843 0.017237 1    -2   1    0   0      #Mock 1
#./HOD_NFW_V14_more_BVG 11.069 0.01213 0.021702 1     0  0.25  0   0      #Mock 2
#./HOD_NFW_V14_more_BVG 11.759 0.08799 0.057523 1.5   0   1    0   0      #Mock 3
#./HOD_NFW_V14_more_BVG 11.759 0.08799 0.057523 1     0   1   500 200     #Mock 4
#./HOD_NFW_V14_more_BVG 11.240 0.02003 0.029333 1     0  0.4   0   0      #Mock 6
#./HOD_NFW_V14_more_BVG 11.264 0.02153 0.030546 1     0  0.25  0   0      #Mock 13
#./HOD_NFW_V14_more_BVG 10.500 0.00636 0.007098 1     0  0.25  0   0      #Mock 16
#./HOD_NFW_V14_more_BVG 11.069 0.01213 0.021702 1     0  0.25  0   0      #Mock 11 reference

Lbox = float(sys.argv[1]) #Mpc/h       #If UNITSIM boxes: 1000
Om  = float(sys.argv[2])               #If UNITSIM cosmology: Om = 0.3089
z1  = float(sys.argv[3])               #If ELG eBOSS: z = 0.865
muAcAs = sys.argv[4]            #Example: mu10.841_Ac0.0062_As0.01418
beta = float(sys.argv[5])              #We consider -2 NI  ; -1-0 B ; 0 P ; 0-... NB

def H(z):
    return H_0*np.sqrt(Om*(1+z)**3+1-Om)

#Galaxy mocks generated 
#bias comprobation of mocks 15/11/2021
data = np.loadtxt('../DESI_outputs/output_V1/mocks_chisquared/galaxies_1000Mpc_V1.4more_NFW_%s_vfact1.00_beta%.3f_K1.00_vt0pm0_BVG_product.dat' % (muAcAs,beta))

X = data[:,0]
Y = data[:,1]
Z = data[:,2]
VZ = data[:,5]
Z_RSD = Z + (1+z1)*VZ/H(z1)

for i in range(0,len(Z_RSD)):
    if Z_RSD[i] < 0:
        Z_RSD[i] = Z_RSD[i]+Lbox
    if Z_RSD[i] > Lbox:
        Z_RSD[i] = Z_RSD[i]-Lbox

out = np.array([X,Y,Z,Z_RSD]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/Mock_%s_beta%.3f.txt' % (muAcAs,beta),out)


#Figura 7 Avila 2020
'''
NB_g = np.loadtxt('../DESI_outputs/output_V1/pruebas/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0089_As0.01954_vfact1.00_beta0.100_K1.00_vt0pm0_BVG_gamma.dat')

X_NB_g = NB_g[:,0]
Y_NB_g = NB_g[:,1]
Z_NB_g = NB_g[:,2]
VZ_NB_g = NB_g[:,5]
Z_RSD_NB_g = Z_NB_g+(1+z1)*VZ_NB_g/H(z1)

for i in range(0,len(Z_RSD_NB_g)):
    if Z_RSD_NB_g[i] < 0:
        Z_RSD_NB_g[i] = Z_RSD_NB_g[i]+Lbox
    if Z_RSD_NB_g[i] > Lbox:
        Z_RSD_NB_g[i] = Z_RSD_NB_g[i]-Lbox

out_NB_g = np.array([X_NB_g,Y_NB_g,Z_NB_g,Z_RSD_NB_g]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/NB_g.txt',out_NB_g)


NB_p = np.loadtxt('../DESI_outputs/output_V1/pruebas/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0089_As0.01954_vfact1.00_beta0.100_K1.00_vt0pm0_BVG_product.dat')

X_NB_p = NB_p[:,0]
Y_NB_p = NB_p[:,1]
Z_NB_p = NB_p[:,2]
VZ_NB_p = NB_p[:,5]
Z_RSD_NB_p = Z_NB_p+(1+z1)*VZ_NB_p/H(z1)

for i in range(0,len(Z_RSD_NB_p)):
    if Z_RSD_NB_p[i] < 0:
        Z_RSD_NB_p[i] = Z_RSD_NB_p[i]+Lbox
    if Z_RSD_NB_p[i] > Lbox:
        Z_RSD_NB_p[i] = Z_RSD_NB_p[i]-Lbox

out_NB_p = np.array([X_NB_p,Y_NB_p,Z_NB_p,Z_RSD_NB_p]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/NB_p.txt',out_NB_p)
'''

'''
m1 = np.loadtxt('../DESI_outputs/output_V1/mocks_fig7_Av2020/galaxies_1000Mpc_V1.4more_NFW_mu10.944_Ac0.0084_As0.01724_vfact1.00_beta-2.000_K1.00_vt0pm0_BVG.dat')

X_m1 = m1[:,0]
Y_m1 = m1[:,1]
Z_m1 = m1[:,2]
VZ_m1 = m1[:,5]
Z_RSD_m1 = Z_m1+(1+z1)*VZ_m1/H(z1)

for i in range(0,len(Z_RSD_m1)):
    if Z_RSD_m1[i] < 0:
        Z_RSD_m1[i] = Z_RSD_m1[i]+Lbox
    if Z_RSD_m1[i] > Lbox:
        Z_RSD_m1[i] = Z_RSD_m1[i]-Lbox

out_m1 = np.array([X_m1,Y_m1,Z_m1,Z_RSD_m1]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/mock1.txt',out_m1)


m2 = np.loadtxt('../DESI_outputs/output_V1/mocks_fig7_Av2020/galaxies_1000Mpc_V1.4more_NFW_mu11.069_Ac0.0121_As0.02170_vfact1.00_beta0.000_K0.25_vt0pm0_BVG.dat')

X_m2 = m2[:,0]
Y_m2 = m2[:,1]
Z_m2 = m2[:,2]
VZ_m2 = m2[:,5]
Z_RSD_m2 = Z_m2+(1+z1)*VZ_m2/H(z1)

for i in range(0,len(Z_RSD_m2)):
    if Z_RSD_m2[i] < 0:
        Z_RSD_m2[i] = Z_RSD_m2[i]+Lbox
    if Z_RSD_m2[i] > Lbox:
        Z_RSD_m2[i] = Z_RSD_m2[i]-Lbox

out_m2 = np.array([X_m2,Y_m2,Z_m2,Z_RSD_m2]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/mock2.txt',out_m2)

m3 = np.loadtxt('../DESI_outputs/output_V1/mocks_fig7_Av2020/galaxies_1000Mpc_V1.4more_NFW_mu11.759_Ac0.0880_As0.05752_vfact1.50_beta0.000_K1.00_vt0pm0_BVG.dat')

X_m3 = m3[:,0]
Y_m3 = m3[:,1]
Z_m3 = m3[:,2]
VZ_m3 = m3[:,5]
Z_RSD_m3 = Z_m3+(1+z1)*VZ_m3/H(z1)

for i in range(0,len(Z_RSD_m3)):
    if Z_RSD_m3[i] < 0:
        Z_RSD_m3[i] = Z_RSD_m3[i]+Lbox
    if Z_RSD_m3[i] > Lbox:
        Z_RSD_m3[i] = Z_RSD_m3[i]-Lbox

out_m3 = np.array([X_m3,Y_m3,Z_m3,Z_RSD_m3]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/mock3.txt',out_m3)

m4 = np.loadtxt('../DESI_outputs/output_V1/mocks_fig7_Av2020/galaxies_1000Mpc_V1.4more_NFW_mu11.759_Ac0.0880_As0.05752_vfact1.00_beta0.000_K1.00_vt500pm200_BVG.dat')

X_m4 = m4[:,0]
Y_m4 = m4[:,1]
Z_m4 = m4[:,2]
VZ_m4 = m4[:,5]
Z_RSD_m4 = Z_m4+(1+z1)*VZ_m4/H(z1)

for i in range(0,len(Z_RSD_m4)):
    if Z_RSD_m4[i] < 0:
        Z_RSD_m4[i] = Z_RSD_m4[i]+Lbox
    if Z_RSD_m4[i] > Lbox:
        Z_RSD_m4[i] = Z_RSD_m4[i]-Lbox

out_m4 = np.array([X_m4,Y_m4,Z_m4,Z_RSD_m4]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/mock4.txt',out_m4)

m5 = np.loadtxt('../DESI_outputs/output_V1/mocks_fig7_Av2020/galaxies_1000Mpc_V1.4more_NFW_mu11.240_Ac0.0200_As0.02933_vfact1.00_beta0.000_K0.40_vt0pm0_BVG.dat')

X_m5 = m5[:,0]
Y_m5 = m5[:,1]
Z_m5 = m5[:,2]
VZ_m5 = m5[:,5]
Z_RSD_m5 = Z_m5+(1+z1)*VZ_m5/H(z1)

for i in range(0,len(Z_RSD_m5)):
    if Z_RSD_m5[i] < 0:
        Z_RSD_m5[i] = Z_RSD_m5[i]+Lbox
    if Z_RSD_m5[i] > Lbox:
        Z_RSD_m5[i] = Z_RSD_m5[i]-Lbox

out_m5 = np.array([X_m5,Y_m5,Z_m5,Z_RSD_m5]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/mock6.txt',out_m5)

m6 = np.loadtxt('../DESI_outputs/output_V1/mocks_fig7_Av2020/galaxies_1000Mpc_V1.4more_NFW_mu11.264_Ac0.0215_As0.03055_vfact1.00_beta0.000_K0.25_vt0pm0_BVG.dat')

X_m6 = m6[:,0]
Y_m6 = m6[:,1]
Z_m6 = m6[:,2]
VZ_m6 = m6[:,5]
Z_RSD_m6 = Z_m6+(1+z1)*VZ_m6/H(z1)

for i in range(0,len(Z_RSD_m6)):
    if Z_RSD_m6[i] < 0:
        Z_RSD_m6[i] = Z_RSD_m6[i]+Lbox
    if Z_RSD_m6[i] > Lbox:
        Z_RSD_m6[i] = Z_RSD_m6[i]-Lbox

out_m6 = np.array([X_m6,Y_m6,Z_m6,Z_RSD_m6]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/mock13.txt',out_m6)

m7 = np.loadtxt('../DESI_outputs/output_V1/mocks_fig7_Av2020/galaxies_1000Mpc_V1.4more_NFW_mu10.500_Ac0.0064_As0.00710_vfact1.00_beta0.000_K0.25_vt0pm0_BVG.dat')

X_m7 = m7[:,0]
Y_m7 = m7[:,1]
Z_m7 = m7[:,2]
VZ_m7 = m7[:,5]
Z_RSD_m7 = Z_m7+(1+z1)*VZ_m7/H(z1)

for i in range(0,len(Z_RSD_m7)):
    if Z_RSD_m7[i] < 0:
        Z_RSD_m7[i] = Z_RSD_m7[i]+Lbox
    if Z_RSD_m7[i] > Lbox:
        Z_RSD_m7[i] = Z_RSD_m7[i]-Lbox

out_m7 = np.array([X_m7,Y_m7,Z_m7,Z_RSD_m7]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/mock16.txt',out_m7)
'''
'''
m8 = np.loadtxt('../DESI_outputs/output_V1/mocks_fig7_Av2020/galaxies_1000Mpc_V1.4more_NFW_mu11.759_Ac0.0880_As0.05752_vfact1.00_beta0.000_K1.00_vt500pm200_BVG.dat')

X_m8 = m8[:,0]
Y_m8 = m8[:,1]
Z_m8 = m8[:,2]
VZ_m8 = m8[:,5]
Z_RSD_m8 = Z_m8+(1+z1)*VZ_m8/H(z1)

for i in range(0,len(Z_RSD_m8)):
    if Z_RSD_m8[i] < 0:
        Z_RSD_m8[i] = Z_RSD_m8[i]+Lbox
    if Z_RSD_m8[i] > Lbox:
        Z_RSD_m8[i] = Z_RSD_m8[i]-Lbox

out_m8 = np.array([X_m8,Y_m8,Z_m8,Z_RSD_m8]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/mock11.txt',out_m8)
'''
'''
m9 = np.loadtxt('../DESI_outputs/output_V1/mocks_fig7_Av2020/galaxies_1000Mpc_V1.4more_NFW_mu11.143_Ac0.0151_As0.02478_vfact1.00_beta0.100_K0.15_vt0pm0_BVG.dat')

X_m9 = m9[:,0]
Y_m9 = m9[:,1]
Z_m9 = m9[:,2]
VZ_m9 = m9[:,5]
Z_RSD_m9 = Z_m9+(1+z1)*VZ_m9/H(z1)

for i in range(0,len(Z_RSD_m9)):
    if Z_RSD_m9[i] < 0:
        Z_RSD_m9[i] = Z_RSD_m9[i]+Lbox
    if Z_RSD_m9[i] > Lbox:
        Z_RSD_m9[i] = Z_RSD_m9[i]-Lbox

out_m9 = np.array([X_m9,Y_m9,Z_m9,Z_RSD_m9]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/mock9.txt',out_m9)
'''

#Comprobacion a partir de mock1 de Avila, 2020 con NI,P,B,NB variando unicamente el parametro beta

'''
NI = np.loadtxt('../DESI_outputs/output_V1/mocks_chisquared/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0089_As0.01954_vfact1.00_beta-2.000_K1.00_vt0pm0_BVG.dat')

X_NI = NI[:,0]
Y_NI = NI[:,1]
Z_NI = NI[:,2]
VZ_NI = NI[:,5]
Z_RSD_NI = Z_NI+(1+z1)*VZ_NI/H(z1)

for i in range(0,len(Z_RSD_NI)):
    if Z_RSD_NI[i] < 0:
        Z_RSD_NI[i] = Z_RSD_NI[i]+Lbox
    if Z_RSD_NI[i] > Lbox:
        Z_RSD_NI[i] = Z_RSD_NI[i]-Lbox

out_NI = np.array([X_NI,Y_NI,Z_NI,Z_RSD_NI]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/n_25_beta-2.0.txt',out_NI)


B_1 = np.loadtxt('../DESI_outputs/output_V1/mocks_chisquared/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0089_As0.01954_vfact1.00_beta-1.000_K1.00_vt0pm0_BVG.dat')

X_B_1 = B_1[:,0]
Y_B_1 = B_1[:,1]
Z_B_1 = B_1[:,2]
VZ_B_1 = B_1[:,5]
Z_RSD_B_1 = Z_B_1+(1+z1)*VZ_B_1/H(z1)

for i in range(0,len(Z_RSD_B_1)):
    if Z_RSD_B_1[i] < 0:
        Z_RSD_B_1[i] = Z_RSD_B_1[i]+Lbox
    if Z_RSD_B_1[i] > Lbox:
        Z_RSD_B_1[i] = Z_RSD_B_1[i]-Lbox

out_B_1 = np.array([X_B_1,Y_B_1,Z_B_1,Z_RSD_B_1]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/n_25_beta-1.0.txt',out_B_1)

B_05 = np.loadtxt('../DESI_outputs/output_V1/mocks_chisquared/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0089_As0.01954_vfact1.00_beta-0.500_K1.00_vt0pm0_BVG.dat')

X_B_05 = B_05[:,0]
Y_B_05 = B_05[:,1]
Z_B_05 = B_05[:,2]
VZ_B_05 = B_05[:,5]
Z_RSD_B_05 = Z_B_05+(1+z1)*VZ_B_05/H(z1)

for i in range(0,len(Z_RSD_B_05)):
    if Z_RSD_B_05[i] < 0:
        Z_RSD_B_05[i] = Z_RSD_B_05[i]+Lbox
    if Z_RSD_B_05[i] > Lbox:
        Z_RSD_B_05[i] = Z_RSD_B_05[i]-Lbox

out_B_05 = np.array([X_B_05,Y_B_05,Z_B_05,Z_RSD_B_05]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/n_25_beta-0.5.txt',out_B_05)

B_025 = np.loadtxt('../DESI_outputs/output_V1/mocks_chisquared/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0089_As0.01954_vfact1.00_beta-0.250_K1.00_vt0pm0_BVG.dat')

X_B_025 = B_025[:,0]
Y_B_025 = B_025[:,1]
Z_B_025 = B_025[:,2]
VZ_B_025 = B_025[:,5]
Z_RSD_B_025 = Z_B_025+(1+z1)*VZ_B_025/H(z1)

for i in range(0,len(Z_RSD_B_025)):
    if Z_RSD_B_025[i] < 0:
        Z_RSD_B_025[i] = Z_RSD_B_025[i]+Lbox
    if Z_RSD_B_025[i] > Lbox:
        Z_RSD_B_025[i] = Z_RSD_B_025[i]-Lbox

out_B_025 = np.array([X_B_025,Y_B_025,Z_B_025,Z_RSD_B_025]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/n_25_beta-0.25.txt',out_B_025)

B_02 = np.loadtxt('../DESI_outputs/output_V1/mocks_chisquared/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0089_As0.01954_vfact1.00_beta-0.200_K1.00_vt0pm0_BVG.dat')

X_B_02 = B_02[:,0]
Y_B_02 = B_02[:,1]
Z_B_02 = B_02[:,2]
VZ_B_02 = B_02[:,5]
Z_RSD_B_02 = Z_B_02+(1+z1)*VZ_B_02/H(z1)

for i in range(0,len(Z_RSD_B_02)):
    if Z_RSD_B_02[i] < 0:
        Z_RSD_B_02[i] = Z_RSD_B_02[i]+Lbox
    if Z_RSD_B_02[i] > Lbox:
        Z_RSD_B_02[i] = Z_RSD_B_02[i]-Lbox

out_B_02 = np.array([X_B_02,Y_B_02,Z_B_02,Z_RSD_B_02]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/n_25_beta-0.2.txt',out_B_02)

B_01 = np.loadtxt('../DESI_outputs/output_V1/mocks_chisquared/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0089_As0.01954_vfact1.00_beta-0.100_K1.00_vt0pm0_BVG.dat')

X_B_01 = B_01[:,0]
Y_B_01 = B_01[:,1]
Z_B_01 = B_01[:,2]
VZ_B_01 = B_01[:,5]
Z_RSD_B_01 = Z_B_01+(1+z1)*VZ_B_01/H(z1)

for i in range(0,len(Z_RSD_B_01)):
    if Z_RSD_B_01[i] < 0:
        Z_RSD_B_01[i] = Z_RSD_B_01[i]+Lbox
    if Z_RSD_B_01[i] > Lbox:
        Z_RSD_B_01[i] = Z_RSD_B_01[i]-Lbox

out_B_01 = np.array([X_B_01,Y_B_01,Z_B_01,Z_RSD_B_01]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/n_25_beta-0.1.txt',out_B_01)

P = np.loadtxt('../DESI_outputs/output_V1/mocks_chisquared/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0089_As0.01954_vfact1.00_beta0.000_K1.00_vt0pm0_BVG.dat')

X_P = P[:,0]
Y_P = P[:,1]
Z_P = P[:,2]
VZ_P = P[:,5]
Z_RSD_P = Z_P+(1+z1)*VZ_P/H(z1)

for i in range(0,len(Z_RSD_P)):
    if Z_RSD_P[i] < 0:
        Z_RSD_P[i] = Z_RSD_P[i]+Lbox
    if Z_RSD_P[i] > Lbox:
        Z_RSD_P[i] = Z_RSD_P[i]-Lbox

out_P = np.array([X_P,Y_P,Z_P,Z_RSD_P]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/n_25_beta0.0.txt',out_P)

NB_01 = np.loadtxt('../DESI_outputs/output_V1/mocks_chisquared/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0089_As0.01954_vfact1.00_beta0.100_K1.00_vt0pm0_BVG.dat')

X_NB_01 = NB_01[:,0]
Y_NB_01 = NB_01[:,1]
Z_NB_01 = NB_01[:,2]
VZ_NB_01 = NB_01[:,5]
Z_RSD_NB_01 = Z_NB_01+(1+z1)*VZ_NB_01/H(z1)

for i in range(0,len(Z_RSD_NB_01)):
    if Z_RSD_NB_01[i] < 0:
        Z_RSD_NB_01[i] = Z_RSD_NB_01[i]+Lbox
    if Z_RSD_NB_01[i] > Lbox:
        Z_RSD_NB_01[i] = Z_RSD_NB_01[i]-Lbox

out_NB_01 = np.array([X_NB_01,Y_NB_01,Z_NB_01,Z_RSD_NB_01]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/n_25_beta0.1.txt',out_NB_01)

NB_02 = np.loadtxt('../DESI_outputs/output_V1/mocks_chisquared/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0089_As0.01954_vfact1.00_beta0.200_K1.00_vt0pm0_BVG.dat')

X_NB_02 = NB_02[:,0]
Y_NB_02 = NB_02[:,1]
Z_NB_02 = NB_02[:,2]
VZ_NB_02 = NB_02[:,5]
Z_RSD_NB_02 = Z_NB_02+(1+z1)*VZ_NB_02/H(z1)

for i in range(0,len(Z_RSD_NB_02)):
    if Z_RSD_NB_02[i] < 0:
        Z_RSD_NB_02[i] = Z_RSD_NB_02[i]+Lbox
    if Z_RSD_NB_02[i] > Lbox:
        Z_RSD_NB_02[i] = Z_RSD_NB_02[i]-Lbox

out_NB_02 = np.array([X_NB_02,Y_NB_02,Z_NB_02,Z_RSD_NB_02]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/n_25_beta0.2.txt',out_NB_02)

NB_03 = np.loadtxt('../DESI_outputs/output_V1/mocks_chisquared/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0089_As0.01954_vfact1.00_beta0.300_K1.00_vt0pm0_BVG.dat')

X_NB_03 = NB_03[:,0]
Y_NB_03 = NB_03[:,1]
Z_NB_03 = NB_03[:,2]
VZ_NB_03 = NB_03[:,5]
Z_RSD_NB_03 = Z_NB_03+(1+z1)*VZ_NB_03/H(z1)

for i in range(0,len(Z_RSD_NB_03)):
    if Z_RSD_NB_03[i] < 0:
        Z_RSD_NB_03[i] = Z_RSD_NB_03[i]+Lbox
    if Z_RSD_NB_03[i] > Lbox:
        Z_RSD_NB_03[i] = Z_RSD_NB_03[i]-Lbox

out_NB_03 = np.array([X_NB_03,Y_NB_03,Z_NB_03,Z_RSD_NB_03]).T
np.savetxt('../DESI_outputs/NERSC_tables/october2021/n_25_beta0.3.txt',out_NB_03)
'''


'''
#COMPROBATION MOCKS WITH BINOMIAL JUNE 2020

n_25_mock1_P = np.loadtxt('../DESI_outputs/output_V1/n_25e-4/galaxies_1000Mpc_V1.4_NFW_mu10.954_Ac0.0089_As0.01954_vfact1.00_beta0.000_K1.00_vt0pm0_2beta0.000_mock_newNB.dat')

X_25_P = n_25_mock1_P[:,0]
Y_25_P = n_25_mock1_P[:,1]
Z_25_P = n_25_mock1_P[:,2]
VZ_25_P = n_25_mock1_P[:,5]
Z_RSD_25_P = Z_25_P+(1+z1)*VZ_25_P/H(z1)

for i in range(0,len(Z_RSD_25_P)):
    if Z_RSD_25_P[i] < 0:
        Z_RSD_25_P[i] = Z_RSD_25_P[i]+Lbox
    if Z_RSD_25_P[i] > Lbox:
        Z_RSD_25_P[i] = Z_RSD_25_P[i]-Lbox

out_25_P = np.array([X_25_P,Y_25_P,Z_25_P,Z_RSD_25_P]).T
np.savetxt('../DESI_outputs/NERSC_tables/n_25_mock1_Poisson.txt',out_25_P)


n_25_mock1_NB = np.loadtxt('../DESI_outputs/output_V1/n_25e-4/galaxies_1000Mpc_V1.4_NFW_mu10.954_Ac0.0089_As0.01954_vfact1.00_beta0.100_K1.00_vt0pm0_2beta0.000_mock_newNB.dat')

X_25_NB = n_25_mock1_NB[:,0]
Y_25_NB = n_25_mock1_NB[:,1]
Z_25_NB = n_25_mock1_NB[:,2]
VZ_25_NB = n_25_mock1_NB[:,5]
Z_RSD_25_NB = Z_25_NB+(1+z1)*VZ_25_NB/H(z1)

for i in range(0,len(Z_RSD_25_NB)):
    if Z_RSD_25_NB[i] < 0:
        Z_RSD_25_NB[i] = Z_RSD_25_NB[i]+Lbox
    if Z_RSD_25_NB[i] > Lbox:
        Z_RSD_25_NB[i] = Z_RSD_25_NB[i]-Lbox

out_25_NB = np.array([X_25_NB,Y_25_NB,Z_25_NB,Z_RSD_25_NB]).T
np.savetxt('../DESI_outputs/NERSC_tables/n_25_mock1_Negbinomial.txt',out_25_NB)



n_25_mock1_B = np.loadtxt('../DESI_outputs/output_V1/n_25e-4/galaxies_1000Mpc_V1.4_NFW_mu10.954_Ac0.0089_As0.01954_vfact1.00_beta0.000_K1.00_vt0pm0_2beta0.100_mock_newNB.dat')

X_25_B = n_25_mock1_B[:,0]
Y_25_B = n_25_mock1_B[:,1]
Z_25_B = n_25_mock1_B[:,2]
VZ_25_B = n_25_mock1_B[:,5]
Z_RSD_25_B = Z_25_B+(1+z1)*VZ_25_B/H(z1)

for i in range(0,len(Z_RSD_25_B)):
    if Z_RSD_25_B[i] < 0:
        Z_RSD_25_B[i] = Z_RSD_25_B[i]+Lbox
    if Z_RSD_25_B[i] > Lbox:
        Z_RSD_25_B[i] = Z_RSD_25_B[i]-Lbox

out_25_B = np.array([X_25_B,Y_25_B,Z_25_B,Z_RSD_25_B]).T
np.savetxt('../DESI_outputs/NERSC_tables/n_25_mock1_binomial.txt',out_25_B)


n_25_mock1_NI = np.loadtxt('../DESI_outputs/output_V1/n_25e-4/galaxies_1000Mpc_V1.4_NFW_mu10.954_Ac0.0089_As0.01954_vfact1.00_beta-2.000_K1.00_vt0pm0_2beta0.000_mock_newNB.dat')

X_25_NI = n_25_mock1_NI[:,0]
Y_25_NI = n_25_mock1_NI[:,1]
Z_25_NI = n_25_mock1_NI[:,2]
VZ_25_NI = n_25_mock1_NI[:,5]
Z_RSD_25_NI = Z_25_NI+(1+z1)*VZ_25_NI/H(z1)

for i in range(0,len(Z_RSD_25_NI)):
    if Z_RSD_25_NI[i] < 0:
        Z_RSD_25_NI[i] = Z_RSD_25_NI[i]+Lbox
    if Z_RSD_25_NI[i] > Lbox:
        Z_RSD_25_NI[i] = Z_RSD_25_NI[i]-Lbox

out_25_NI = np.array([X_25_NI,Y_25_NI,Z_25_NI,Z_RSD_25_NI]).T
np.savetxt('../DESI_outputs/NERSC_tables/n_25_mock1_Nearestinteger.txt',out_25_NI)



'''
#MOCKS DECEMBER 2020
'''
n_25_mock1 = np.loadtxt('../DESI_outputs/output_V1/n_25e-4/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0089_As0.01954_vfact1.00_beta-2.000_K1.00_vt0pm0_mock.dat')

X_25_1 = n_25_mock1[:,0]
Y_25_1 = n_25_mock1[:,1]
Z_25_1 = n_25_mock1[:,2]
VZ_25_1 = n_25_mock1[:,5]
Z_RSD_25_1 = Z_25_1+(1+z1)*VZ_25_1/H(z1)

for i in range(0,len(Z_RSD_25_1)):
    if Z_RSD_25_1[i] < 0:
        Z_RSD_25_1[i] = Z_RSD_25_1[i]+Lbox
    if Z_RSD_25_1[i] > Lbox:
        Z_RSD_25_1[i] = Z_RSD_25_1[i]-Lbox

out_25_1 = np.array([X_25_1,Y_25_1,Z_25_1,Z_RSD_25_1]).T
np.savetxt('../DESI_outputs/NERSC_tables/n_25_mock1_mod.txt',out_25_1)



n_25_mock2 = np.loadtxt('../DESI_outputs/output_V1/n_25e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.079_Ac0.0129_As0.02471_vfact1.00_beta0.000_K0.25_vt0pm0_mock.dat')

X_25_2 = n_25_mock2[:,0]
Y_25_2 = n_25_mock2[:,1]
Z_25_2 = n_25_mock2[:,2]
VZ_25_2 = n_25_mock2[:,5]
Z_RSD_25_2 = Z_25_2+(1+z1)*VZ_25_2/H(z1)

for i in range(0,len(Z_RSD_25_2)):
    if Z_RSD_25_2[i] < 0:
        Z_RSD_25_2[i] = Z_RSD_25_2[i]+Lbox
    if Z_RSD_25_2[i] > Lbox:
        Z_RSD_25_2[i] = Z_RSD_25_2[i]-Lbox

out_25_2 = np.array([X_25_2,Y_25_2,Z_25_2,Z_RSD_25_2]).T
np.savetxt('../DESI_outputs/NERSC_tables/n_25_mock2_mod.txt',out_25_2)



n_25_mock3 = np.loadtxt('../DESI_outputs/output_V1/n_25e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.723_Ac0.0860_As0.05916_vfact1.50_beta0.000_K1.00_vt0pm0_mock.dat')

X_25_3 = n_25_mock3[:,0]
Y_25_3 = n_25_mock3[:,1]
Z_25_3 = n_25_mock3[:,2]
VZ_25_3 = n_25_mock3[:,5]
Z_RSD_25_3 = Z_25_3+(1+z1)*VZ_25_3/H(z1)

for i in range(0,len(Z_RSD_25_3)):
    if Z_RSD_25_3[i] < 0:
        Z_RSD_25_3[i] = Z_RSD_25_3[i]+Lbox
    if Z_RSD_25_3[i] > Lbox:
        Z_RSD_25_3[i] = Z_RSD_25_3[i]-Lbox

out_25_3 = np.array([X_25_3,Y_25_3,Z_25_3,Z_RSD_25_3]).T
np.savetxt('../DESI_outputs/NERSC_tables/n_25_mock3_mod.txt',out_25_3)



n_25_mock4 = np.loadtxt('../DESI_outputs/output_V1/n_25e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.723_Ac0.0860_As0.05916_vfact1.00_beta0.000_K1.00_vt500pm200_mock.dat')

X_25_4 = n_25_mock4[:,0]
Y_25_4 = n_25_mock4[:,1]
Z_25_4 = n_25_mock4[:,2]
VZ_25_4 = n_25_mock4[:,5]
Z_RSD_25_4 = Z_25_4+(1+z1)*VZ_25_4/H(z1)

for i in range(0,len(Z_RSD_25_4)):
    if Z_RSD_25_4[i] < 0:
        Z_RSD_25_4[i] = Z_RSD_25_4[i]+Lbox
    if Z_RSD_25_4[i] > Lbox:
        Z_RSD_25_4[i] = Z_RSD_25_4[i]-Lbox

out_25_4 = np.array([X_25_4,Y_25_4,Z_25_4,Z_RSD_25_4]).T
np.savetxt('../DESI_outputs/NERSC_tables/n_25_mock4_mod.txt',out_25_4)



n_25_mock9 = np.loadtxt('../DESI_outputs/output_V1/n_25e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.153_Ac0.0160_As0.02829_vfact1.00_beta0.100_K0.15_vt0pm0_mock.dat')

X_25_5 = n_25_mock9[:,0]
Y_25_5 = n_25_mock9[:,1]
Z_25_5 = n_25_mock9[:,2]
VZ_25_5 = n_25_mock9[:,5]
Z_RSD_25_5 = Z_25_5+(1+z1)*VZ_25_5/H(z1)

for i in range(0,len(Z_RSD_25_5)):
    if Z_RSD_25_5[i] < 0:
        Z_RSD_25_5[i] = Z_RSD_25_5[i]+Lbox
    if Z_RSD_25_5[i] > Lbox:
        Z_RSD_25_5[i] = Z_RSD_25_5[i]-Lbox

out_25_5 = np.array([X_25_5,Y_25_5,Z_25_5,Z_RSD_25_5]).T
np.savetxt('../DESI_outputs/NERSC_tables/n_25_mock9_mod.txt',out_25_5)




















n_20_mock1 = np.loadtxt('../DESI_outputs/output_V1/n_20e-4/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0071_As0.01563_vfact1.00_beta-2.000_K1.00_vt0pm0_mock.dat')

X_20_1 = n_20_mock1[:,0]
Y_20_1 = n_20_mock1[:,1]
Z_20_1 = n_20_mock1[:,2]
VZ_20_1 = n_20_mock1[:,5]
Z_RSD_20_1 = Z_20_1+(1+z1)*VZ_20_1/H(z1)

for i in range(0,len(Z_RSD_20_1)):
    if Z_RSD_20_1[i] < 0:
        Z_RSD_20_1[i] = Z_RSD_20_1[i]+Lbox
    if Z_RSD_20_1[i] > Lbox:
        Z_RSD_20_1[i] = Z_RSD_20_1[i]-Lbox

out_20_1 = np.array([X_20_1,Y_20_1,Z_20_1,Z_RSD_20_1]).T
np.savetxt('../DESI_outputs/NERSC_tables/n_20_mock1_mod.txt',out_20_1)



n_20_mock2 = np.loadtxt('../DESI_outputs/output_V1/n_20e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.079_Ac0.0103_As0.01977_vfact1.00_beta0.000_K0.25_vt0pm0_mock.dat')

X_20_2 = n_20_mock2[:,0]
Y_20_2 = n_20_mock2[:,1]
Z_20_2 = n_20_mock2[:,2]
VZ_20_2 = n_20_mock2[:,5]
Z_RSD_20_2 = Z_20_2+(1+z1)*VZ_20_2/H(z1)

for i in range(0,len(Z_RSD_20_2)):
    if Z_RSD_20_2[i] < 0:
        Z_RSD_20_2[i] = Z_RSD_20_2[i]+Lbox
    if Z_RSD_20_2[i] > Lbox:
        Z_RSD_20_2[i] = Z_RSD_20_2[i]-Lbox

out_20_2 = np.array([X_20_2,Y_20_2,Z_20_2,Z_RSD_20_2]).T
np.savetxt('../DESI_outputs/NERSC_tables/n_20_mock2_mod.txt',out_20_2)



n_20_mock3 = np.loadtxt('../DESI_outputs/output_V1/n_20e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.723_Ac0.0688_As0.04733_vfact1.50_beta0.000_K1.00_vt0pm0_mock.dat')

X_20_3 = n_20_mock3[:,0]
Y_20_3 = n_20_mock3[:,1]
Z_20_3 = n_20_mock3[:,2]
VZ_20_3 = n_20_mock3[:,5]
Z_RSD_20_3 = Z_20_3+(1+z1)*VZ_20_3/H(z1)

for i in range(0,len(Z_RSD_20_3)):
    if Z_RSD_20_3[i] < 0:
        Z_RSD_20_3[i] = Z_RSD_20_3[i]+Lbox
    if Z_RSD_20_3[i] > Lbox:
        Z_RSD_20_3[i] = Z_RSD_20_3[i]-Lbox

out_20_3 = np.array([X_20_3,Y_20_3,Z_20_3,Z_RSD_20_3]).T
np.savetxt('../DESI_outputs/NERSC_tables/n_20_mock3_mod.txt',out_20_3)




n_20_mock4 = np.loadtxt('../DESI_outputs/output_V1/n_20e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.723_Ac0.0688_As0.04733_vfact1.00_beta0.000_K1.00_vt500pm200_mock.dat')

X_20_4 = n_20_mock4[:,0]
Y_20_4 = n_20_mock4[:,1]
Z_20_4 = n_20_mock4[:,2]
VZ_20_4 = n_20_mock4[:,5]
Z_RSD_20_4 = Z_20_4+(1+z1)*VZ_20_4/H(z1)

for i in range(0,len(Z_RSD_20_4)):
    if Z_RSD_20_4[i] < 0:
        Z_RSD_20_4[i] = Z_RSD_20_4[i]+Lbox
    if Z_RSD_20_4[i] > Lbox:
        Z_RSD_20_4[i] = Z_RSD_20_4[i]-Lbox

out_20_4 = np.array([X_20_4,Y_20_4,Z_20_4,Z_RSD_20_4]).T
np.savetxt('../DESI_outputs/NERSC_tables/n_20_mock4_mod.txt',out_20_4)



n_20_mock9 = np.loadtxt('../DESI_outputs/output_V1/n_20e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.153_Ac0.0128_As0.02263_vfact1.00_beta0.100_K0.15_vt0pm0_mock.dat')

X_20_5 = n_20_mock9[:,0]
Y_20_5 = n_20_mock9[:,1]
Z_20_5 = n_20_mock9[:,2]
VZ_20_5 = n_20_mock9[:,5]
Z_RSD_20_5 = Z_20_5+(1+z1)*VZ_20_5/H(z1)

for i in range(0,len(Z_RSD_20_5)):
    if Z_RSD_20_5[i] < 0:
        Z_RSD_20_5[i] = Z_RSD_20_5[i]+Lbox
    if Z_RSD_20_5[i] > Lbox:
        Z_RSD_20_5[i] = Z_RSD_20_5[i]-Lbox

out_20_5 = np.array([X_20_5,Y_20_5,Z_20_5,Z_RSD_20_5]).T
np.savetxt('../DESI_outputs/NERSC_tables/n_20_mock9_mod.txt',out_20_5)





n_5_mock1 = np.loadtxt('../DESI_outputs/output_V1/n_5e-4/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0018_As0.00391_vfact1.00_beta-2.000_K1.00_vt0pm0_mock.dat')

X_5_1 = n_5_mock1[:,0]
Y_5_1 = n_5_mock1[:,1]
Z_5_1 = n_5_mock1[:,2]
VZ_5_1 = n_5_mock1[:,5]
Z_RSD_5_1 = Z_5_1+(1+z1)*VZ_5_1/H(z1)

for i in range(0,len(Z_RSD_5_1)):
    if Z_RSD_5_1[i] < 0:
        Z_RSD_5_1[i] = Z_RSD_5_1[i]+Lbox
    if Z_RSD_5_1[i] > Lbox:
        Z_RSD_5_1[i] = Z_RSD_5_1[i]-Lbox

out_5_1 = np.array([X_5_1,Y_5_1,Z_5_1,Z_RSD_5_1]).T
np.savetxt('../DESI_outputs/NERSC_tables/n_5_mock1_mod.txt',out_5_1)



n_5_mock2 = np.loadtxt('../DESI_outputs/output_V1/n_5e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.079_Ac0.0026_As0.00494_vfact1.00_beta0.000_K0.25_vt0pm0_mock.dat')

X_5_2 = n_5_mock2[:,0]
Y_5_2 = n_5_mock2[:,1]
Z_5_2 = n_5_mock2[:,2]
VZ_5_2 = n_5_mock2[:,5]
Z_RSD_5_2 = Z_5_2+(1+z1)*VZ_5_2/H(z1)

for i in range(0,len(Z_RSD_5_2)):
    if Z_RSD_5_2[i] < 0:
        Z_RSD_5_2[i] = Z_RSD_5_2[i]+Lbox
    if Z_RSD_5_2[i] > Lbox:
        Z_RSD_5_2[i] = Z_RSD_5_2[i]-Lbox

out_5_2 = np.array([X_5_2,Y_5_2,Z_5_2,Z_RSD_5_2]).T
np.savetxt('../DESI_outputs/NERSC_tables/n_5_mock2_mod.txt',out_5_2)



n_5_mock3 = np.loadtxt('../DESI_outputs/output_V1/n_5e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.723_Ac0.0172_As0.01183_vfact1.50_beta0.000_K1.00_vt0pm0_mock.dat')

X_5_3 = n_5_mock3[:,0]
Y_5_3 = n_5_mock3[:,1]
Z_5_3 = n_5_mock3[:,2]
VZ_5_3 = n_5_mock3[:,5]
Z_RSD_5_3 = Z_5_3+(1+z1)*VZ_5_3/H(z1)

for i in range(0,len(Z_RSD_5_3)):
    if Z_RSD_5_3[i] < 0:
        Z_RSD_5_3[i] = Z_RSD_5_3[i]+Lbox
    if Z_RSD_5_3[i] > Lbox:
        Z_RSD_5_3[i] = Z_RSD_5_3[i]-Lbox

out_5_3 = np.array([X_5_3,Y_5_3,Z_5_3,Z_RSD_5_3]).T
np.savetxt('../DESI_outputs/NERSC_tables/n_5_mock3_mod.txt',out_5_3)



n_5_mock4 = np.loadtxt('../DESI_outputs/output_V1/n_5e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.723_Ac0.0172_As0.01183_vfact1.00_beta0.000_K1.00_vt500pm200_mock.dat')

X_5_4 = n_5_mock4[:,0]
Y_5_4 = n_5_mock4[:,1]
Z_5_4 = n_5_mock4[:,2]
VZ_5_4 = n_5_mock4[:,5]
Z_RSD_5_4 = Z_5_4+(1+z1)*VZ_5_4/H(z1)

for i in range(0,len(Z_RSD_5_4)):
    if Z_RSD_5_4[i] < 0:
        Z_RSD_5_4[i] = Z_RSD_5_4[i]+Lbox
    if Z_RSD_5_4[i] > Lbox:
        Z_RSD_5_4[i] = Z_RSD_5_4[i]-Lbox

out_5_4 = np.array([X_5_4,Y_5_4,Z_5_4,Z_RSD_5_4]).T
np.savetxt('../DESI_outputs/NERSC_tables/n_5_mock4_mod.txt',out_5_4)



n_5_mock9 = np.loadtxt('../DESI_outputs/output_V1/n_5e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.153_Ac0.0032_As0.00566_vfact1.00_beta0.100_K0.15_vt0pm0_mock.dat')

X_5_5 = n_5_mock9[:,0]
Y_5_5 = n_5_mock9[:,1]
Z_5_5 = n_5_mock9[:,2]
VZ_5_5 = n_5_mock9[:,5]
Z_RSD_5_5 = Z_5_5+(1+z1)*VZ_5_5/H(z1)

for i in range(0,len(Z_RSD_5_5)):
    if Z_RSD_5_5[i] < 0:
        Z_RSD_5_5[i] = Z_RSD_5_5[i]+Lbox
    if Z_RSD_5_5[i] > Lbox:
        Z_RSD_5_5[i] = Z_RSD_5_5[i]-Lbox

out_5_5 = np.array([X_5_5,Y_5_5,Z_5_5,Z_RSD_5_5]).T
np.savetxt('../DESI_outputs/NERSC_tables/n_5_mock9_mod.txt',out_5_5)
'''

