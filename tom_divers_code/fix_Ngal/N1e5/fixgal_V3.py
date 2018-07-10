import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



def f(N, Ngal_):
	return N-Ngal


DLM = np.array([0.00, 1.00, 2.00])
M1 = np.array([10.00, 10.50, 11.00, 11.50, 12.00, 12.50, 13.00, 13.50, 14.00])
Ngal= 100000
tol=Ngal*1e-2



for i in range(0, len(DLM)):
	for j in range(0, len(M1)):
		print("*",i,j,"*")
		Ll=12.0
		Lh= 16.0
		Lmid=(Ll+Lh)/2
		FILENAME= "../../HOD_Santi/output/galaxies_V1.0_M1_%.2f_DLM%.2f_Lth%.2f.dat" % (M1[j], DLM[i], Ll)
		array= np.loadtxt(FILENAME, delimiter=' ')
		N= len(array[:, 7])
		print("Nfound=",N)
		Nmid=N
		if ( (N<= Ngal- tol) ):
			print('Not enough Galaxies')
		else: 			
			while abs(f(Nmid, Ngal)) > tol:
				Nmid=np.size(np.where(array[:,7]>Lmid))
				Nl=np.size(np.where(array[:,7]>Ll))
				Nh=np.size(np.where(array[:,7]>Lh))
				print("     N(",Ll,")=",Nl)
				print("     N(",Lmid,")=",Nmid)
				print("     N(",Lh,")=",Nh)
				if (f(Nl, Ngal)*f(Nmid, Ngal)>0):
					Ll = Lmid
				else:
					Lh = Lmid
				Lmid=(Ll+Lh)/2
				Nmid=np.size(np.where(array[:,7]>Lmid))

			print(Lmid)
			out = array[np.where(array[:,7]>Lmid)]

			np.savetxt("galaxies_V1.0_M1_%.2f_DLM%.2f.dat"%(M1[j], DLM[i]),out)
			 

