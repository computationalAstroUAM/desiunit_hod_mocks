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

#1.    THEORETICAL PK
#definimos la cosmologia considerada en UNITSIM:
Om = 0.3089
Omb = 0.049
h = 0.6774
ns = 0.9667

#elegimos las escalas en k. No hace falta elegir mucho rango. Estamos asumiendo bias lineal, a grandes escalas no lo es.
Ngrid = 512
BoxSize = 1000.

Fundamental =  2*np.pi/BoxSize
knyquist = np.pi*Ngrid/BoxSize

kmin = Fundamental
kmax = knyquist
dk = 0.007
k = np.arange(kmin,kmax,dk)

#1.1.   Using camb, definimos los parametros cosmologicos
params = cb.set_params(ns=ns,H0=h*100,ombh2=Omb*h**2,omch2=(Om-Omb)*h**2,WantTransfer=True) #Set all CAMB parameters at once, including parameters which are part of the CAMBparams structure, as well as global parameters.
PK_th = get_matter_power_interpolator(params,kmax = kmax,nonlinear=False,hubble_units=True)

#z1 = 0.9873
z1 = 0.8594
Pk_th = PK_th.P(z1,k)

#normalizacion
sigma8 = 0.8147
def sigma2_R(P,k,R):
    return spi.simps(P*1./(2*pi**2)*k**2*((3./(k*R))*spherical_jn(1,k*R))**2,x=k)

def Omega_M(z):
    return (Om*(1+z)**3)/(Om*(1+z)**3+1-Om)

gamma = 0.545 #LCDM + GR
def f(z):
    return (Omega_M(z))**gamma

def integrand(z):
    return f(z)/(1+z)

def D(z): #Growth factor
    return np.exp(integrate.fixed_quad(integrand, z,0)[0]) #el [1] de esta funcion te da el error


sigma8_lin = np.sqrt(sigma2_R(Pk_th,k,8))
Pk_th = (sigma8**2/sigma8_lin**2)*Pk_th * D(z1)**2
#print((sigma8**2/sigma8_lin**2)*D(z1)**2)
#print(D(z1)**2)
#print(sigma8_lin)



#1.2.    Using UNIT particles 

catalogue = '../../DESI_outputs/dm_particles_0.005_100'
names = ['x','y','z']
f = CSVCatalog(catalogue, names)
numhalos = len(f['x'])
f['Position'] = f['x'][:, None]*[1,0,0]+ f['y'][:, None]*[0,1,0]+ f['z'][:, None]*[0,0,1]

f.attrs['BoxSize'] = BoxSize #anadimos a f el tamano de la caja

mesh = f.to_mesh (compensated = True, Nmesh = Ngrid , BoxSize = BoxSize , position = "Position" )
r = FFTPower(mesh, mode='1d',dk=dk,kmin=Fundamental-dk/2) # se calcula el power spectrum en una dimension con la FFT. En mesh le anadimos los datos que nenesita de la linea anterior
Pk = r.power
Pk_th_particles = Pk['power'].real-Pk.attrs['shotnoise']









#2.   PK UNITSIM haloes
catalogue = sys.argv[1] 


#Verbose = False   #Verbose is a general programming term for produce lots of logging output. You can think of it as asking the program to tell me everything about what you are doing all the time. Just set it to true and see what happens.


#print('Nyquist frequency:',knyquist)
names = ['x','y','z']
f = CSVCatalog(catalogue, names) 
numhalos = len(f['x'])
f['Position'] = f['x'][:, None]*[1,0,0]+ f['y'][:, None]*[0,1,0]+ f['z'][:, None]*[0,0,1]

f.attrs['BoxSize'] = BoxSize #anadimos a f el tamano de la caja

mesh = f.to_mesh (compensated = True, Nmesh = Ngrid , BoxSize = BoxSize , position = "Position" )
r = FFTPower(mesh, mode='1d',dk=dk,kmin=Fundamental-dk/2) # se calcula el power spectrum en una dimension con la FFT. En mesh le anadimos los datos que nenesita de la linea anterior
Pk = r.power
#print("Shotnoise",Pk.attrs['shotnoise'], "=?",BoxSize**3/numhalos)
#Shotnoise = BoxSize**3/numhalos
#print(Shotnoise)
Pk_sim = Pk['power'].real-Pk.attrs['shotnoise']

n = numhalos/BoxSize**3

sell = np.where((k < 0.1) & (k > 0.01))
#sell = np.where((k < 0.1) & (k > 0.015))

k = k[sell]
Pk_sim = Pk_sim[sell]
Pk_th = Pk_th[sell]
PK_th_particles = Pk_th_particles[sell]
#print(k)
#print(Pk['k'][sell])
#print(Pk_sim)

DeltaP2 = ((2*np.pi)**2/(k**2*dk*Vbox))*(Pk_sim+1/n)**2

b_min = 0
b_max = 12
db = 0.002
num = (b_max-b_min)/db
num = int(num)

xi_2 = np.zeros((num))

b_min_range = int(b_min/db)
b_max_range = int(b_max/db)
i=0
for i in range(b_min_range,b_max_range):
    for j in range(0,len(Pk_th)):
        xi_2[i] += ((((i*db)**2.)*Pk_th_particles[j]-Pk_sim[j])**2.)/DeltaP2[j]
    i+=1

bias = np.arange(b_min,b_max+0.001,db)

pos = np.where(xi_2 == np.amin(xi_2))
Bias = bias[pos]
Bias = float(Bias)
print(Bias)

#np.savetxt('bias/bias_75.txt',Bias)

#fig, axs = plt.subplots(2, 1,sharex=True)
#fig.subplots_adjust(hspace=0)
#axs[0].loglog(k, Pk_sim,'b-',label=r'$P_{UNIT}$')
#axs[0].loglog(k, Pk_th*Bias**2,'g-',label=r'$P_{th}\cdot b^2$')
#axs[0].fill_between(k, Pk_sim - np.sqrt(DeltaP2), Pk_sim + np.sqrt(DeltaP2),color='yellow', alpha=0.3)
#axs[0].set_xlabel(r'$k$')
#axs[0].set_ylabel(r'$P(k)$')
#axs[0].text(0.02,900,r'$14 < log_{10}\left(M\right) < 14.2$')
#axs[0].text(0.02,2000,r'$z = 0.9873$')
#axs[0].legend()

#axs[1].plot(k,Pk_sim/(Pk_th*Bias**2),'r-',label=r'$\frac{P_{UNIT}}{P_{th}\cdot b^2}$')
#axs[1].set_xlabel(r'$k$')
#axs[1].fill_between(k, Pk_sim/(Pk_th*Bias**2) - np.sqrt(DeltaP2)/(Pk_th*Bias**2), Pk_sim/(Pk_th*Bias**2) + np.sqrt(DeltaP2)/(Pk_th*Bias**2),color='yellow', alpha=0.4)
#axs[1].axhline(y=1,color='k')
#axs[1].legend()
#plt.savefig('plots/Pk_75_bin_k_01.png',dpi=300)
#plt.close()



#Bias calculation  with 2PCF comparing with particles   November 2021


particles = np.loadtxt('../../FCFC/isotropic_2PCF_redshiftspace/xi_snap100_bin5.txt') 

r = particles[:,0]
xi0 = particles[:,3]

#haloes = #2PCF from haloes snapshot 100




