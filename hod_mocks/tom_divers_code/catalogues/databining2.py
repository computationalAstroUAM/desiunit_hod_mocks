import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
matplotlib.rcParams.update({'font.size': 12})

FILENAME= "OuterRim_STEP266_500_XYZM_cent.dat"
array= np.loadtxt(FILENAME)
logM=np.log10(array[:, 3])
#print(array)
delM= 0.5
Mh = np.arange(14, 15, delM)

for j in range(0, len(Mh)):
	i = np.where(np.logical_and(logM > Mh[j], logM<(Mh[j]+delM)))[0]
#	print(i)
#	print(Mh[j], Mh[j]+delM)
	newarray = array[i, :]
	print(newarray)
	np.savetxt("4biasplots/OuterRim_STEP266_500_cXYZM_{}logMh{}.dat".format(Mh[j], Mh[j]+delM), newarray)

