import numpy as np
import sys
import io


#edges = np.concatenate( [ np.array( np.arange(10.58,12.4999,0.04)), np.array(np.arange(12.5,13.6,0.05)),  np.array(np.arange(13.6,14.0,0.1)),np.array(np.arange(14.2,14.6,0.2))])

#We have UNITSIM catalogues created expressely with 21 particles which corresponds to 10**10.418 


edges = np.concatenate( [ np.array( np.arange(10.4,12.4,0.04)), np.array(np.arange(12.4,13.6,0.05)),  np.array(np.arange(13.6,14.0,0.1)),np.array(np.arange(14.2,14.6,0.2)),np.array([15.00001])]) 


X,Y,Z,M,ID = np.loadtxt('../../DESI_outputs/UNIT_001/halos_UNIT_001/outp_100_4096_x_y_z_vx_vy_vz_m_id_no_subhalos_logM_mayor_10.418.txt',usecols = (0,1,2,6,7),unpack = True) # El catalogo debe tener datos de X, Y, Z. Velocidades no es necesario al trabajar en real space

#For z=0.9873:  out_97p_X_Y_Z_VX_VY_VZ_M.txt
#For z=0.8594:  out_100p_X_Y_Z_VX_VY_VZ_M.txt

#sel0 = np.where((ID == -1))

#M = M[sel0]
#X = X[sel0]
#Y = Y[sel0]
#Z = Z[sel0]

for i in range(0,len(edges)-1):
    sel = np.where((M < edges[i+1]) & (M > edges[i]))
    pos = np.array([X[sel],Y[sel],Z[sel]]).T
    nombre = str(i) + ".txt"
    np.savetxt('../../DESI_outputs/positions_for_each_mass_bin_PNG_snap100_21particles/' + nombre,pos)


