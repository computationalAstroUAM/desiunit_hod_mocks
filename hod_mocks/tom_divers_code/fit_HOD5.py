import numpy as np
from scipy.special import erf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import matplotlib.colors as colors

def delta_sq(b_mod, b_targ, nmod, n_targ):
        return (((b_mod**2 - b_targ**2)/btarg**2 )**2 + ((nmod - n_targ)/n_targ)**2)

def Trap(fa,fb,h):
        return h/2 * (fa+fb)

def CalcNsat(M, M0, M1, alpha):
	if (M*alpha>M0*alpha):
        	return (((10**M) - 10**M0)/10**M1)**alpha
	else:
		0.

def CalcNcent(M, Mmin, sig):
	return 0.5*( 1 + erf((M - Mmin)/sig))

def my_arange(a, b, dr, decimals=6):
	res = [a]
	k=1
	while res[-1] < b:
		tmp = round(a+k*dr, decimals)
		if tmp >b:
			break
		res.append(tmp)
		k+=1
	return np.asarray(res)

M0 = my_arange(11.0 ,15.5 ,0.02, 2)
Mmin = my_arange(11.0 ,15.5 , 0.02, 2)
alpha  = -0.5
sig = 1.0
V= 500.0**3
ngal= np.zeros((len(M0),len(Mmin)))
bgal= np.zeros((len(M0),len(Mmin)))
delt = np.zeros_like(bgal)

#Favole
btarg = 1.25
ntarg = 4.0e-4

#eboss
#btarg = 1.25
#ntarg = 0.000283385216095

print "## Traget ##"
print "b=",btarg, " n=",ntarg
print " "

File = "Scripts_Tom/CUTE_box/plots/Halos/Mh_n_b_Halostest_halfres.txt"
array = np.loadtxt(File)
n = array[:, 1]
bh = array[:, 2]
logMh = array[:, 0]
Dsave=np.zeros((len(M0)*len(Mmin)))
nsave=np.zeros((len(M0)*len(Mmin)))
bsave=np.zeros((len(M0)*len(Mmin)))
M0save=np.zeros((len(M0)*len(Mmin)))
Mminsave=np.zeros((len(M0)*len(Mmin)))
count = 0

for k in range(0,len(M0)):
	for l in range(0,len(Mmin)):
		if (M0[k]>Mmin[l]):
			M1 = Mmin[l]+1.0
			I_N= np.zeros_like(logMh)
			I_b= np.zeros_like(logMh)
			I_ngal= np.zeros_like(logMh)
			#initialize
			Nsat = np.zeros_like(logMh)
			Ncent= np.zeros_like(logMh)
			if logMh[0]>M0[k]:
				Nsat[0] = CalcNsat(logMh[0], M0[k], M1, alpha)
			Ncent[0] = CalcNcent(logMh[0],Mmin[l], sig)
			fN= np.zeros_like(logMh)
			fb=np.zeros_like(logMh)
			fN[0]=(Nsat[0]+Ncent[0])
			fb[0] = fN[0] * bh[0]
			for i in range(1,len(logMh)):
				delM= logMh[i]- logMh[i-1]
				if logMh[i]>M0[k]:
					Nsat[i] = CalcNsat(logMh[i], M0[k], M1, alpha)
				Ncent[i]= CalcNcent(logMh[i], Mmin[l], sig)
				fN[i] = (Nsat[i]+Ncent[i])
				I_N[i-1]= Trap(fN[i-1],fN[i], delM)*n[i-1]
				I_b[i-1]= Trap(fN[i-1],fN[i], delM)*n[i-1]*bh[i-1]
				I_ngal[i-1]= I_N[i-1]/(V )
			Ngal = np.sum(I_N)
			bgal[k, l] = np.sum(I_b)/Ngal
			ngal[k,l] = np.sum(I_ngal)
			delt[k,l]= delta_sq(bgal[k,l], btarg, ngal[k,l], ntarg)
                        #if (delt[k,l]<0.85
			if ((bgal[k, l]<1.5*btarg) & (bgal[k, l]>0.5*btarg) & (ngal[k,l]<10*ntarg) & (ngal[k,l]>ntarg*0.1)):
                                print("Delta=",delt[k,l]," b=",bgal[k,l]," n=",ngal[k,l],"M0=",M0[k],"Mmin=",Mmin[l])
                                bsave[count]=bgal[k,l]
                                nsave[count]=ngal[k,l]
                                M0save[count]=M0[k]
                                Mminsave[count]=Mmin[l]
                                Dsave[count]=delt[k,l]
                                count = count +1
                        #if (delt[k,l]<0.8):
                        #iif M0[k]==10.24 and Mmin[l]==10.67:
                        #       print(ngal[k,l], bgal[k,l])
		else:
			delt[k,l]=10.

#print(bgal)

#biasplo
#for k in range(0, len(M0)):
#	for l in range(0, len(Mmin)):
#		if (M0[k]>Mmin[l]) and bgal[k,l]>1.2 :
#			print(M0[k],Mmin[l], bgal[k,l])
#			plt.scatter(Mmin[l], M0[k], c=(bgal[k,l]), s=60)#, cmap='viridis', s= 60)# ,norm=colors.LogNorm(vmin=np.min(bgal), vmax=np.max(bgal)))
			
#plt.colorbar()
m,n= np.where(delt==np.min(delt))
#,n= np.where(delt==np.minimum(delt))
#print(np.min(np.absolute(delt)))
#print('bgal=',bgal[m,n])
#print('diff=', abs(bgal[m,n]-1.1))
#plt.plot(Mmin[m][0],M0[n][0], 'k*', markersize=20, alpha=0.5)
#plt.xlabel('Mmin')
#plt.xlim(14.6 ,15.0)
#plt.ylabel('M0')
#plt.ylim(14.8 ,15.2)
#plt.title('Galaxy Bias')
#plt.savefig('parameter_bias22.png')

print(bgal[-1,-1], 'last')



#number density
#for k in range(0, len(M0)):
#	for l in range(0, len(Mmin)):
#		if (M0[k]>Mmin[l]):
#			plt.scatter(Mmin[l], M0[k], c=(ngal[k,l]) , cmap=cm.jet, s=60, vmin=np.min(ngal), vmax=np.max(ngal))
#



#print('Mimimum:')
#k=m
#l=n
#print("Delta=",delt[k,l][0]," b=",bgal[k,l][0]," n=",ngal[k,l][0],"M0=",M0[k][0],"Mmin=",Mmin[l][0])


#print('diff=', ngal[m,n]-1.1)
#print('M0=',M0[m], 'Mmin=',Mmin[n])
#plt.plot(Mmin[l][0],M0[k][0], 'k*', markersize=20, alpha=0.5)
#plt.xlabel('Mmin')
#plt.xlim(14.7 ,14.9)
#plt.ylabel('M0')
#plt.ylim(14.92 ,15.1)
#plt.title('Galaxy Number Density')
#plt.colorbar()
#plt.savefig('parameter_ndens_3.png')

#print('saves ',count,'entries')

#out = np.transpose(np.array([Dsave[0:count],bsave[0:count],nsave[0:count],M0save[0:count],Mminsave[0:count]]))
#np.savetxt("save.txt", out, header=str("Delta b n M0 Mmin"))

#for k in range(0, len(M0)):
#        for l in range(0, len(Mmin)):
#                if (M0[k]>Mmin[l]):
#                        plt.scatter(Mmin[l], M0[k], c=(bgal[k,l]) , cmap=cm.jet, s=60, vmin=np.min(bgal), vmax=np.max(bgal))
#



print('Mimimum:')
k=m
l=n
print("Delta=",delt[k,l][0]," b=",bgal[k,l][0]," n=",ngal[k,l][0],"M0=",M0[k][0],"Mmin=",Mmin[l][0])


#print('diff=', ngal[m,n]-1.1)
#print('M0=',M0[m], 'Mmin=',Mmin[n])
#plt.plot(Mmin[l][0],M0[k][0], 'k*', markersize=20, alpha=0.5)
#plt.xlabel('Mmin')
#plt.xlim(14.7 ,14.9)
#plt.ylabel('M0')
#plt.ylim(14.92 ,15.1)
#plt.title('Galaxy Bias')
#plt.colorbar()
#plt.savefig('parameter_BIAS.png')
                                                              
