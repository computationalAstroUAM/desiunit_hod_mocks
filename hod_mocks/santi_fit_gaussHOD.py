import numpy as np
from scipy.special import erf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import matplotlib.colors as colors

#Metrics that we want to minimisize (or "error")
def delta_sq(b_mod, b_targ, nmod, n_targ):
        return (((b_mod**2 - b_targ**2)/btarg**2 )**2 + ((nmod - n_targ)/n_targ)**2)

def Trap(fa,fb,h):
        return h/2 * (fa+fb)

#<Nsat(Mh)>
def CalcNsat(M, M0, M1, alpha,fsat):
        if M>M0:
                return fsat*(((10**M) - 10**M0)/10**M1)**alpha
        else:
                return 0.

#<Ncen(Mh)>
def CalcNcent(M, mu, sig,fcen):
	return fcen/(sig*np.sqrt(np.pi*2.0))*np.exp(-(M-mu)**2/(2.0*sig**2))
	#return fcen/(sig*np.sqrt(np.pi*2.0))*np.exp(-(10**M-10**mu)**2/(2.0*sig**2))



#parameters we loop through (fitting grid)
f = np.arange(0.002, 0.1 ,0.002) 
Mmean = np.arange(11.0,15.0,0.01)


#scaling params
fsat = f/6.0
M0 = Mmean - 0.1
M1 = Mmean + 0.3
#fixed params
alpha  = 0.8
sig = 0.12 


#Volume which was used to computed N(Mh)  (=n(Mh) * V)
V= 500.0**3

ngal= np.zeros((len(f),len(Mmean)))
bgal= np.zeros((len(f),len(Mmean)))
delt = np.zeros_like(bgal)


#eboss ELGV2
ntarg = 0.000246356087788
btarg = 1.22

print "## Target ##"
print "b=",btarg, " n=",ntarg
print " "

#Mh, N(Mh), b(Mh)
File = "../CUTE_box/plots/Halos/Mh_n_b_Halostest_halfres.txt"
array = np.loadtxt(File)
nh = array[:, 1]
bh = array[:, 2]
logMh = array[:, 0]

Dsave=np.zeros((len(f)*len(Mmean)))
nsave=np.zeros((len(f)*len(Mmean)))
bsave=np.zeros((len(f)*len(Mmean)))
fsave=np.zeros((len(f)*len(Mmean)))
Mmeansave=np.zeros((len(f)*len(Mmean)))
count = 0


#loop through the paramters (Acen and mu)
#and compute (n,b) 
for k in range(0,len(f)):
	for l in range(0,len(Mmean)):
		I_N= np.zeros_like(logMh)
		I_b= np.zeros_like(logMh)
		I_ngal= np.zeros_like(logMh)
		#initialize
		Nsat = np.zeros_like(logMh)
		Ncent= np.zeros_like(logMh)
		if (logMh[0]>M0[l] and fsat[k]!=0.):
			Nsat[0] = CalcNsat(logMh[0], M0[l], M1[l], alpha,fsat[k])
		Ncent[0] = CalcNcent(logMh[0],Mmean[l], sig,f[k])
		fN= np.zeros_like(logMh)
		fb=np.zeros_like(logMh)
		fN[0]=(Nsat[0]+Ncent[0])
		fb[0] = fN[0] * bh[0]
		for i in range(1,len(logMh)):
			delM= logMh[i]- logMh[i-1]
			if (logMh[i]>M0[l] and fsat[k]!=0.):
				Nsat[i] = CalcNsat(logMh[i], M0[l], M1[l], alpha,fsat[k])
			Ncent[i]= CalcNcent(logMh[i], Mmean[l], sig,f[k])
			fN[i] = (Nsat[i]+Ncent[i])
			I_N[i-1]= Trap(fN[i-1],fN[i], delM)*nh[i-1]
			#I_Nsat[i-1]= Trap(fN[i-1],Nsat[i], delM)*nh[i-1]
			I_b[i-1]= Trap(fN[i-1],fN[i], delM)*nh[i-1]*bh[i-1]
			I_ngal[i-1]= I_N[i-1]/(V )
		Ngal = np.sum(I_N)
		bgal[k, l] = np.sum(I_b)/Ngal
		ngal[k,l] = np.sum(I_ngal)
		#compute the "error" wrt the target quantities
		delt[k,l]= delta_sq(bgal[k,l], btarg, ngal[k,l], ntarg)


#find the parameters that give us the minimum "error"
print "delt=",delt
m,n= np.where(delt==np.min(delt))
print('Mimimum:')
k=m
l=n
#print("Delta=",delt[k,l])
print("Delta=",delt[k,l][0]," b=",bgal[k,l][0]," n=",ngal[k,l][0],"f=",f[k][0],"Mmean=",Mmean[l][0])
print("Fix params: ","alpha=",alpha," sigma=",sig)
print("re-scaled params: fsat=",f[k][0]/6.0,"M0=",Mmean[l][0] - 0.1,"M1=",Mmean[l][0] + 0.3)


#print(bgal)

#biasplo
#for k in range(0, len(M0)):
#	for l in range(0, len(Mmin)):
#		if (M0[k]>Mmin[l]) and bgal[k,l]>1.2 :
#			print(M0[k],Mmin[l], bgal[k,l])
#			plt.scatter(Mmin[l], M0[k], c=(bgal[k,l]), s=60)#, cmap='viridis', s= 60)# ,norm=colors.LogNorm(vmin=np.min(bgal), vmax=np.max(bgal)))
			
#plt.colorbar()

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



