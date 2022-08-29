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


Vbox = 10**9. #Mpc^3/h^3
#z1 = 0.9873
z1 = 0.8594


#1.    THEORETICAL PK
#definimos la cosmologia considerada en UNITSIM:
Om = 0.3089
Omb = 0.049
h = 0.6774
ns = 0.9667

#elegimos las escalas en k. No hace falta elegir mucho rango. Estamos asumiendo bias lineal, a grandes escalas no lo es.
#Ngrid = 512
#BoxSize = 1000.

#Fundamental =  2*np.pi/BoxSize
#knyquist = np.pi*Ngrid/BoxSize
#dk = Fundamental

#kmin = Fundamental+dk/2.
#kmax = knyquist+dk/2.
#k = np.arange(kmin,kmax,dk)

def Omega_M(z):
    return (Om*(1+z)**3)/(Om*(1+z)**3+1-Om)

gamma = 0.545 #LCDM + GR
def f(z):
    return (Omega_M(z))**gamma

#1. Pk using UNIT particles 
k,Pk_th_particles = np.loadtxt('../../DESI_outputs/Pk_particles/Pk_particles.txt',usecols = (0,1),unpack = True) #We used previously dm_particles_0.005_100, dm_particles_0.5_100 dm_particles_0.5_100_InvPhase_tsc.Pk

#if we want to use theory of camb unit:
k,Pk_th_particles = np.loadtxt('../Pk_theory/Pk_camb_unit.txt',usecols = (0,1),unpack = True)


#2.   PK UNITSIM haloes
Pk_number = sys.argv[1]
dat = np.loadtxt(Pk_number)
Pk_sim = dat[:,1]
n = dat[0,3]

sell = np.where((k < 0.1) & (k > 0.01))
#sell = np.where((k < 0.1) & (k > 0.015))


k = k[sell]
Pk_sim = Pk_sim[sell]
Pk_th_particles = Pk_th_particles[sell]

dk = k[1]-k[0]
#print(Pk_th)
#print(Pk_th_particles)
#print(len(Pk_th))
#print(len(Pk_th_particles))
#print(Pk_sim)
#print(len(Pk_sim))

DeltaP2 = ((2*np.pi)**2/(k**2*dk*Vbox))*(Pk_sim+1/n)**2

b_min = 0
b_max = 20
db = 0.002
num = (b_max-b_min)/db
num = int(num)

xi_2 = np.zeros((num))

b_min_range = int(b_min/db)
b_max_range = int(b_max/db)
i=0
for i in range(b_min_range,b_max_range):
    for j in range(0,len(Pk_th_particles)):
        xi_2[i] += ((((i*db)**2.)*Pk_th_particles[j]-Pk_sim[j])**2.)/DeltaP2[j] #real space
        #xi_2[i] += ((((i*db)**2.+ (2/3.)*(i*db)*f(z1) + (1/5.)*f(z1)**2)*Pk_th_particles[j]-Pk_sim[j])**2.)/DeltaP2[j]  #redshift space
    i+=1

bias = np.arange(b_min,b_max+0.001,db)
Bias = float(bias[np.where(xi_2 == np.amin(xi_2))])
print(Bias)


#1. hacer dos subsets de xi_2, uno con los valores despues del minimo y otro con los valores antes del minimo


posicion_min = int(np.where(xi_2==min(xi_2))[0])
xi_2_antes = np.zeros((posicion_min))
#print('posicion_min',posicion_min)

for i in range(0,posicion_min):
    xi_2_antes[i] = xi_2[i]

xi_2_despues = np.zeros((len(xi_2)-posicion_min-1))
for i in range(posicion_min+1,len(xi_2)):
    xi_2_despues[i-posicion_min-1] = xi_2[i]

#2. para cada subset, buscar la posicion del valor de la chi squared el cual se acerque mas a min(xi_2)+1

a = int(np.where(np.abs(xi_2_antes-np.amin(xi_2)-1) == np.amin(np.abs(xi_2_antes-np.amin(xi_2)-1)))[0])
Bias_antes = bias[a]
print(Bias_antes)

b = posicion_min + 1 + int(np.where(np.abs(xi_2_despues-np.amin(xi_2)-1) == np.amin(np.abs(xi_2_despues-np.amin(xi_2)-1)))[0])
Bias_despues = bias[b]
print(Bias_despues)


