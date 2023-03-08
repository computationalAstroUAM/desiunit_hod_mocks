import numpy as np


bias = np.zeros(82)
bdown = np.zeros(82)
bup = np.zeros(82)
for i in range(0,82):
    bias[i] = np.loadtxt('../../DESI_outputs/bias_PNG_snap100_21particles/bias_%d.txt' % i)[0]
    bdown[i] = np.loadtxt('../../DESI_outputs/bias_PNG_snap100_21particles/bias_%d.txt' % i)[1]
    bup[i] = np.loadtxt('../../DESI_outputs/bias_PNG_snap100_21particles/bias_%d.txt' % i)[2]

print(bias)
print(bdown)
print(bup)

data  = np.loadtxt('../../DESI_outputs/output/HMF_and_accumulated_and_log_mainhalos_snap100_extended_21_particles.txt')
mhalo = data[:,0]
print(mhalo)
print(len(mhalo))
print(len(bias))
out = np.array([mhalo,bias,bdown,bup]).T

np.savetxt('../../DESI_outputs/bias_PNG_snap100_21particles/all_bias.txt',out)


