import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



def f(N, Ngal_):
	return N-Ngal


DLM = np.array([0.00])
M1 = np.array([10.00, 11.00, 12.00, 13.00, 14.00])

n=7.642602e-04
L=500.0
Ngal= int(n * (L*L*L))
tol=Ngal*1e-2



for i in range(0, len(DLM)):
	for j in range(0, len(M1)):
		print("*",i,j,"*")
		Ll=10.0
		Lh= 18.0
		Lmid=(Ll+Lh)/2
		FILENAME= "../../HOD_Santi/output500/galaxies_V1.0_M1_%.2f_DLM%.2f_Lth%.2f.dat" % (M1[j], DLM[i], Ll)
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

			np.savetxt("Santi_HOD_500/galaxiesNL3_V1.0_M1_%.2f_DLM%.2f.dat"%(M1[j], DLM[i]),out)
			 

