import numpy as np
from scipy.special import spherical_jn #importamos la funcion esferica de bessel que necesitaremos despues
import scipy.integrate as spi
import scipy.integrate as integrate


# 1. Analysing effects of the linear approach made in slide 159 Santi slides)

data_mock1 = np.loadtxt('hod_mocks/output_V1/n_25e-4/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0089_As0.01954_vfact1.00_beta-2.000_K1.00_vt0pm0_mock.dat')

VZ_1 = data_mock1[:,5]

vz_max_1 = np.max(VZ_1)

for i in range(0,len(VZ_1)):
    if VZ_1[i]<0:
        VZ_1[i] = (-1.)*VZ_1[i]

vz_mean_1 = np.mean(VZ_1)

z = 0.9873
c = 299792.458 #Km/s

#Relativistic formula that relates redhsift with Hubble recesion velocity
v_Hubble = (((z+1)**2-1)*c)/(1+(z+1)**2)

print('v_Hubble=',v_Hubble)
print('For a first mock we consider maximum radial velocity a mean radial velocity and we compare it to the recession velocity correpondent to z=0.9873')
print('vz_max=',vz_max_1)
print('vz_mean=',vz_mean_1)
print('vz_max/v_Hubble=',vz_max_1/v_Hubble)
print('vz_mean/v_Hubble=',vz_mean_1/v_Hubble)

data_mock2 = np.loadtxt('hod_mocks/output_V1/n_25e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.079_Ac0.0129_As0.02471_vfact1.00_beta0.000_K0.25_vt0pm0_mock.dat')

VZ_2 = data_mock2[:,5]

vz_max_2 = np.max(VZ_2)

for i in range(0,len(VZ_2)):
    if VZ_2[i]<0:
        VZ_2[i] = (-1.)*VZ_2[i]

vz_mean_2 = np.mean(VZ_2)

print('For a second mock we have')
print('vz_max=',vz_max_2)
print('vz_mean=',vz_mean_2)
print('vz_max/v_Hubble=',vz_max_2/v_Hubble)
print('vz_mean/v_Hubble=',vz_mean_2/v_Hubble)



# 2. Analysing the linear growth rate f approximation

a1 = 1/(1+z)
Om = 0.3089
OL = 1+Om
def Omega_M(a):
    return (Om*a**(-3))/(Om*a**(-3)+1-Om)

def Omega_L(a):
    return (OL)/(Om*a**(-3)+OL)


gamma = 0.545 #LCDM + GR

#We define the approximation for f
def f(a):
    return (Omega_M(a))**gamma

print('Approximation for f=',f(a1))

H_0 = 100 #Mpc/h
def H(a):
    return np.sqrt(Om*a**(-3)+1-Om)


from scipy.integrate import quad

def integrand(a):
    return 1./((a**3.)*(H(a))**3.)


I = quad(integrand, a1, 1)[0]


def g(a):
    return (5/2.)*Omega_M(a)*I


f1 = -1-Omega_M(a1)/2.+Omega_L(a1)+1./I
print('More realistic f equation 4 paper Hamilton=',f1)
print('Error in f',f1/f(a1)-1)

#Lahav et al 1991

f2 = Omega_M(a1)**(4/7.)+(1+Omega_M(a1)/2.)*(Omega_L(a1)/70.)

print('Aproximation de Lahav 1991de f=',f2)






