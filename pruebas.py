import sys, os, glob, getpass
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np

T_b_dir = '../'

densityparameter = np.genfromtxt(T_b_dir+'meanbrightnesstemperatureanddensityparameter.txt')


Vbox=10**(9.) #Mpc/h
edges = np.concatenate( [ np.array( np.arange(10.58,12.4999,0.04)), np.array(np.arange(12.5,13.6,0.05)),  np.array(np.arange(13.6,14.0,0.1)),np.array(np.arange(14.2,14.6,0.2))])

dm = edges[1:]-edges[:-1]
mhist = edges[1:]-0.5*dm

M=10**np.linspace(10.4,15,1000)

i=0
hmf = np.zeros(len(mhist),dtype=float)
deltaM = np.zeros(len(mhist),dtype=float)
sel = np.zeros(len(mhist),dtype=float)

for Mi in mhist:
    sel=np.where((edges[i] < np.log10(M)) & (np.log10(M) < edges[i+1]))
    deltaM=10**(edges[i+1])-10**(edges[i])
    hmf[i]=len(M[sel])/(deltaM*Vbox)
    i+=1


#out = np.array([hmf,mhist]).T
#np.savetxt('prueba.txt',out)



'''
print(len(edges))

M = 10**np.linspace(10.5,14.5,2000) 
X = np.linspace(0,1999,2000)
Y = np.linspace(0,1999,2000)
Z = np.linspace(0,1999,2000)
edges = np.array(edges)
#Escribir unos arrays mas grandes donde quepa todo
for i in range(0,4):
    sel = np.where((np.log10(M) < edges[i+1]) & (np.log10(M) > edges[i]))
    pos = np.array([X[sel],Y[sel],Z[sel]]).T
    nombre = str(i) + ".txt"
    np.savetxt('positions_for_each_mass_bin/' + nombre, pos)
'''


a = np.zeros(10)

b = np.array([2,3,4,5,6,7,8])
c = np.array([3,4,5,6,7,8,9])

for i in range(0,len(a)):
    for j in range(0,len(b)):
        a[i] += (i**2)*(c[j]-b[j])**2
    i+=1

#print(a)


data = np.loadtxt('out_97p_X_Y_Z_VX_VY_VZ_M.txt')

ID = data[:,7]
sel = np.where((ID ==-1))

X=data[:,0]
Y=data[:,1]
Z=data[:,2]
VX=data[:,3]
VY=data[:,4]
VZ=data[:,5]
M = data[:,6]

X = X[sel]
Y = Y[sel]
Z = Z[sel]
VX = VX[sel]
VY = VY[sel]
VZ = VZ[sel]
M = M[sel]
m = np.log10(M)
ID = ID[sel]


out = np.array([X,Y,Z,VX,VY,VZ,m,ID]).T
np.savetxt('out_97p_X_Y_Z_VX_VY_VZ_logM.txt',out)












