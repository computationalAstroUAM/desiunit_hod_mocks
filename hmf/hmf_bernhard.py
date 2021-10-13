import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

Vbox=10**(9.) #Mpc/h
edges = np.concatenate( [ np.array( np.arange(10.58,12.4999,0.04)), np.array(np.arange(12.5,13.6,0.05)),  np.array(np.arange(13.6,14.0,0.1)),np.array(np.arange(14.2,14.6,0.2))]) 
dm = edges[1:]-edges[:-1]
mhist = edges[1:]-0.5*dm

data = np.loadtxt('../../DESI_outputs/out_100p_X_Y_Z_VX_VY_VZ_M.txt') 
halos = data[:,7]
sel0 = np.where((halos == -1))  #nos quedamos solo con los halos y no subhalos
M = data[:,6]
M = M[sel0]

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

out=np.array([mhist,hmf,hmf_acc,hmflog]).T
np.savetxt('../../DESI_outputs/output/HMF_and_accumulated_and_log_mainhalos_snap100.txt',out)

