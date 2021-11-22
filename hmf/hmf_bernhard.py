import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np





Vbox=10**(9.) #(Mpc/h)^3
#edges = np.concatenate( [ np.array( np.arange(10.58,12.4999,0.04)), np.array(np.arange(12.5,13.6,0.05)),  np.array(np.arange(13.6,14.0,0.1)),np.array(np.arange(14.2,14.6,0.2))]) 

edges = np.concatenate( [ np.array( np.arange(9.5,12.5,0.04)), np.array(np.arange(12.5,13.6,0.05)),  np.array(np.arange(13.6,14.0,0.1)),np.array(np.arange(14.2,15,0.2))])

dm = edges[1:]-edges[:-1]
mhist = edges[1:]-0.5*dm

data = np.loadtxt('../../DESI_outputs/out_100p_X_Y_Z_VX_VY_VZ_M.txt') 
halos = data[:,7]
sel0 = np.where((halos == -1))  #nos quedamos solo con los halos y no subhalos
M = data[:,6]
M = M[sel0]

#M = 10**np.array([8.001,8.1001,8.2001,8.3001,8.4001])

#saber cuantos halo hay por bin de masas para generar un array de bines adecuado
bines = 0.02
minimo = 8
maximo = 16
A = np.linspace(minimo,maximo,int((maximo-minimo)/bines)+1)
B = np.zeros(len(A),dtype=float)
Sel = np.zeros(len(A),dtype=float)
i=0 
for i in range(0,len(A)-1):
    Sel = np.where((np.log10(M) > A[i]) & (np.log10(M) < A[i+1]))
    B[i] = len(M[Sel])
    i+=1

out = np.array([A,B]).T
np.savetxt('numberofhalos_per_bin.txt',out)
'''
i=0
hmf = np.zeros(len(mhist),dtype=float)
deltaM = np.zeros(len(mhist),dtype=float)
sel = np.zeros(len(mhist),dtype=float)
hmflog = np.zeros(len(mhist),dtype=float)
deltalogM = np.zeros(len(mhist),dtype=float)



for Mi in mhist:
    sel=np.where((edges[i] < np.log10(M)) & (np.log10(M) < edges[i+1]))
    deltaM=10**(edges[i+1])-10**(edges[i])
    deltalogM = edges[i+1]-edges[i]
    hmf[i] = len(M[sel])/(deltaM*Vbox)
    hmflog[i] = len(M[sel])/(deltalogM*Vbox)
    i+=1

k=0
hmf_acc = np.zeros(len(mhist),dtype=float)
sel3 = np.zeros(len(mhist),dtype=float)

for k in range(0,len(mhist)):
    sel3=np.where((mhist[k] < np.log10(M)))
    hmf_acc[k]=len(M[sel3])/Vbox
    k+=1
'''
#out=np.array([mhist,hmf,hmf_acc,hmflog]).T
#np.savetxt('../../DESI_outputs/output/HMF_and_accumulated_and_log_mainhalos_snap100.txt',out)

