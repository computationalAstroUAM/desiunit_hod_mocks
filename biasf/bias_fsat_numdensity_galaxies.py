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

#THEORETICAL PK
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

#utilizando las funciones de camb, definimos los parametros cosmologicos
params = cb.set_params(ns=ns,H0=h*100,ombh2=Omb*h**2,omch2=(Om-Omb)*h**2,WantTransfer=True) #Set all CAMB parameters at once, including parameters which are part of the CAMBparams structure, as well as global parameters.
PK_th = get_matter_power_interpolator(params,kmax = kmax,nonlinear=False,hubble_units=True)
results=cb.get_results(params)
z1 = 0.9873
Pk_th = PK_th.P(0,k)

sigma8_2 = results.get_sigma8()

#normalizacion
sigma8 = 0.8147
def sigma2_R(P,k,R):
    return spi.simps(P*1./(2*pi**2)*k**2*((3./(k*R))*spherical_jn(1,k*R))**2,x=k)

def Omega_M(z):
    return (Om*(1+z)**3)/(Om*(1+z)**3+1-Om)
gamma = 0.545 #LCDM + GR
def ff(z):
    return (Omega_M(z))**gamma

def integrand(z):
    return ff(z)/(1+z)

def D(z): #Growth factor
    return np.exp(integrate.fixed_quad(integrand, z,0)[0]) #el [1] de esta funcion te da el error


sigma8_lin = np.sqrt(sigma2_R(Pk_th,k,8))
print(sigma8_lin)
print(sigma8_2)
Pk_th = (sigma8**2/sigma8_lin**2)*Pk_th * D(z1)**2









catalogue = sys.argv[1]

names1 = ['x','y','z','vx+Dvx','vy+Dvy','vz+Dvz','M','ID','Dvx','Dvy','Dvz','Dx','Dy','Dz','id']
f = CSVCatalog(catalogue, names1)

numgalaxies = len(f['x'])
print('Number of galaxies:',len(f['x']))

ID = f['ID']
sel = np.where((ID ==-1))
ID = np.array(ID)
fsat = len(ID[sel])/len(ID)
print('fsat:',fsat)


#Verbose = False   #Verbose is a general programming term for produce lots of logging output. You can think of it as asking the program to tell me everything about what you are doing all the time. Just set it to true and see what happens.

knyquist = pi*Ngrid/BoxSize #calculamos la frecuencia de Nyquist
print('Nyquist frequency:',knyquist)

numberdensity = numgalaxies/BoxSize**3
print('number density:',numberdensity)

names = ['x','y','z']
f['Position'] = f['x'][:, None]*[1,0,0]+ f['y'][:, None]*[0,1,0]+ f['z'][:, None]*[0,0,1]
f.attrs['BoxSize'] = BoxSize #anadimos a f el tamano de la caja

mesh = f.to_mesh (compensated = True, Nmesh = Ngrid , BoxSize = BoxSize , position = "Position" )
r = FFTPower(mesh, mode='1d',dk=dk,kmin=Fundamental-dk/2,poles=[0,2,4]) # se calcula el power spectrum en una dimension con la FFT. En mesh le anadimos los datos que nenesita de la linea anterior
Pk = r.power
Pk_poles = r.poles
#print("Shotnoise",Pk.attrs['shotnoise'], "=?",BoxSize**3/numgalaxiesi)
#Shotnoise = BoxSize**3/numgalaxies
#print(Shotnoise)
Pk_sim = Pk['power'].real-Pk.attrs['shotnoise']

n = numgalaxies/BoxSize**3

sell = np.where((k < 0.1) & (k > 0.01 ))

k = k[sell]
Pk_sim = Pk_sim[sell]
Pk_th_realspace = Pk_th[sell]

DeltaP2 = ((2*np.pi)**2/(k**2*dk*Vbox))*(Pk_sim+1/n)**2

b_min = 0
b_max = 3
db = 0.002
num = (b_max-b_min)/db
num = int(num)

xi_2 = np.zeros((num))
i=0


b_min_range = int(b_min/db)
b_max_range = int(b_max/db)

for i in range(b_min_range,b_max_range):
    for j in range(0,len(Pk_th_realspace)):
        xi_2[i] += ((((i*db)**2.)*Pk_th_realspace[j]-Pk_sim[j])**2.)/DeltaP2[j]
    i+=1

bias = np.arange(b_min,b_max+0.001,db)

pos = np.where(xi_2 == np.amin(xi_2))
Bias = bias[pos]
print('Bias:',Bias)

P_0 = Pk_poles['power_0'].real -Pk.attrs['shotnoise']  #aqui usamos todo el rango en k
P_2 = Pk_poles['power_2'].real  #aqui usamos todo el rango en k
P_4 = Pk_poles['power_4'].real #aqui usamos todo el rango en k

k = Pk['k']  #aqui usamos todo el rango en k

out = np.array([k,P_0,P_2,P_4]).T
np.savetxt('data_multipoles/n_25e-4_mock1_fsat0_realspace.txt',out)
















#Redshift space distortions

#H_0 = 100 #Mpc/h
#def H(z):
#    return H_0*np.sqrt(Om*(1+z)**3+1-Om)

#f['z_RSD'] = f['z']+(1+z1)*f['vz+Dvz']/H(z1)

#f['Position_redshift'] = f['x'][:, None]*[1,0,0]+ f['y'][:, None]*[0,1,0]+ f['z_RSD'][:, None]*[0,0,1]

#mesh = f.to_mesh (compensated = True, Nmesh = Ngrid , BoxSize = BoxSize , position = "Position_redshift" )
#r = FFTPower(mesh, mode='1d',dk=dk,kmin=Fundamental-dk/2,poles=[0,2,4]) # se calcula el power spectrum en una dimension con la FFT. En mesh le anadimos los datos que nenesita de la linea anterior
#Pk = r.power

#Pk_sim_rsd = Pk['power'].real-Pk.attrs['shotnoise']

#Pk_poles = r.poles

#P_0 = Pk_poles['power_0'].real -Pk.attrs['shotnoise']
#P_2 = Pk_poles['power_2'].real
#P_4 = Pk_poles['power_4'].real

#modes = Pk_poles['modes']

#k = Pk['k']

#out = np.array([k,P_0,P_2,P_4]).T
#np.savetxt('data_multipoles/n_25e-4_mock1_fsat0.txt',out)
#n = numgalaxies/BoxSize**3

#sell = np.where((k < 0.1) & (k > 0.01 ))

#k = k[sell]
#Pk_sim_rsd = Pk_sim_rsd[sell]
#Pk_th = Pk_th[sell]

#print(Pk_sim)

#DeltaP2_rsd = ((2*np.pi)**2/(k**2*dk*Vbox))*(Pk_sim_rsd+1/n)**2

#b_min_rsd = 0
#b_max_rsd = 3
#db_rsd = 0.002
#num_rsd = (b_max_rsd-b_min_rsd)/db_rsd
#num_rsd = int(num_rsd)

#xi_2_rsd = np.zeros((num_rsd))
#i=0


#b_min_range_rsd = int(b_min_rsd/db_rsd)
#b_max_range_rsd = int(b_max_rsd/db_rsd)

#for i in range(b_min_range_rsd,b_max_range_rsd):
#    for j in range(0,len(Pk_th)):
#        xi_2_rsd[i] += ((((i*db_rsd)**2.+(2/3.)*(i*db_rsd)*(ff(z1))+0.2*(ff(z1))**2)*Pk_th[j]-Pk_sim_rsd[j])**2.)/DeltaP2_rsd[j]
#    i+=1

#print(xi_2)
#bias_rsd = np.arange(b_min_rsd,b_max_rsd+0.001,db_rsd)

#pos_rsd = np.where(xi_2_rsd == np.amin(xi_2_rsd))
#Bias_rsd = bias_rsd[pos_rsd]
#print('Bias_RSD:',Bias_rsd)



#KAISER COMPROBATION

kmax = knyquist
data = np.loadtxt('data_multipoles/n_25e-4_mock1.txt')
k = data[:,0]

params = cb.set_params(ns=ns,H0=h*100,ombh2=Omb*h**2,omch2=(Om-Omb)*h**2,WantTransfer=True) #Set all CAMB parameters at once, including parameters which are part of the CAMBparams structure, as well as global parameters.
PK_th = get_matter_power_interpolator(params,kmax = kmax,nonlinear=False,hubble_units=True)


z1 = 0.9873
Pk_th = PK_th.P(z1,k)

#normalizacion
sigma8 = 0.8147
def sigma2_R(P,k,R):
    return spi.simps(P*1./(2*pi**2)*k**2*((3./(k*R))*spherical_jn(1,k*R))**2,x=k)

def Omega_M(z):
    return (Om*(1+z)**3)/(Om*(1+z)**3+1-Om)
gamma = 0.545 #LCDM + GR
def ff(z):
    return (Omega_M(z))**gamma

def integrand(z):
    return ff(z)/(1+z)

def D(z): #Growth factor
    return np.exp(integrate.fixed_quad(integrand, z,0)[0]) #el [1] de esta funcion te da el error


sigma8_lin = np.sqrt(sigma2_R(Pk_th,k,8))
Pk_th = (sigma8**2/sigma8_lin**2)*Pk_th * D(z1)**2

#RSD
P_0_Kaiser = (Bias**2+(2/3.)*ff(z1)*Bias+(1/5.)*(ff(z1))**2)*Pk_th  #Normalmente va Bias_rsd si esta puesto Bias
P_2_Kaiser = ((4/3.)*Bias*ff(z1)+(4/7.)*(ff(z1))**2)*Pk_th
P_4_Kaiser = (8/35.)*(ff(z1))**2*Pk_th
P_0_Kaiser = np.array(P_0_Kaiser)
P_2_Kaiser = np.array(P_2_Kaiser)
P_4_Kaiser = np.array(P_4_Kaiser)

Kaiser = np.array([k,P_0_Kaiser,P_2_Kaiser,P_4_Kaiser]).T

np.savetxt('Kaiser_multipoles/n_25e-4_mock1_fsat0_KAISER_biasrealspace.txt',Kaiser)


#Real space

#P_0_Kaiser_realspace = (Bias**2)*Pk_th  #Normalmente va Bias si esta puesto Bias_rsd
#P_0_Kaiser_realspace = np.array(P_0_Kaiser_realspace)

#Kaiser_realspace = np.array([k,P_0_Kaiser_realspace]).T

#np.savetxt('Kaiser_multipoles/n_25e-4_mock1_fsat0_KAISER_realspace_biasrsd.txt',Kaiser_realspace)









