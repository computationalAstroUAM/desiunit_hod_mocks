import numpy as np
import camb as cb
from camb import get_matter_power_interpolator
from nbodykit.source.catalog import CSVCatalog
from nbodykit.lab import FFTPower #Fast Fourier Transform, tecnica muy estandar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import io


edges = np.concatenate( [ np.array( np.arange(10.58,12.4999,0.04)), np.array(np.arange(12.5,13.6,0.05)),  np.array(np.arange(13.6,14.0,0.1)),np.array(np.arange(14.2,14.6,0.2))])
print(len(edges))
data = np.loadtxt('../../DESI_outputs/out_100p_X_Y_Z_VX_VY_VZ_M.txt') # El catalogo debe tener datos de X, Y, Z. Velocidades no es necesario al trabajar en real space

#For z=0.9873:  out_97p_X_Y_Z_VX_VY_VZ_M.txt
#For z=0.8594:  out_100p_X_Y_Z_VX_VY_VZ_M.txt


sel0 = np.where((data[:,7]==-1))

X = data[:,0]
Y = data[:,1]
Z = data[:,2]
M = data[:,6]

M = M[sel0]
X = X[sel0]
Y = Y[sel0]
Z = Z[sel0]

for i in range(0,len(edges)-1):
    sel = np.where((np.log10(M) < edges[i+1]) & (np.log10(M) > edges[i]))
    pos = np.array([X[sel],Y[sel],Z[sel]]).T
    nombre = str(i) + ".txt"
    np.savetxt('../../DESI_outputs/positions_for_each_mass_bin_snap100/' + nombre,pos)


