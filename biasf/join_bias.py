import numpy as np


bias = np.zeros(77)
for i in range(0,77):
    bias[i] = np.loadtxt('../../DESI_outputs/bias_snap100/bias_particles_short_%d.txt' % i)


data  = np.loadtxt('../../DESI_outputs/output/HMF_and_accumulated_and_log_mainhalos_snap100.txt')
mhalo = data[:,0] 


out = np.array([mhalo,bias]).T

np.savetxt('../../DESI_outputs/bias_snap100/all_bias_particles_short.txt',out)


