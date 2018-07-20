import scipy.integrate as integrate
import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

def Calcgsat(M, M0, M1, alpha):
    Nsat_Mh = 0. 
    if (M>M0):
        Nsat_Mh = (((10**M) - 10**M0)/10**M1)**alpha
    return Nsat_Mh

gsat = np.vectorize(Calcgsat)

def Calcgcent(M, mu, sig):
    return 1.0/(sig*np.sqrt(np.pi*2.0))*np.exp(-(M - mu)**2/(2.0*sig**2))
gcent = np.vectorize(Calcgcent)

def CalcAsat(nd, vol, fsat, Is):
    return nd*vol*fsat/Is

Asat = np.vectorize(CalcAsat)

def CalcAcen(nd, vol, fsat, Ic):
    return nd*vol*(1.0 - fsat)/Ic

Acen = np.vectorize(CalcAcen)

#Metrics that we want to minimisize (or "error")
def delta_sq(b_mod, b_targ):
    delt =(b_mod**2 - b_targ**2)/b_targ**2 )**2
    return delt

#############################                                      
#                                                                  
#  Input ARGUMENTS                                                 
#                                                                  
narg = len(sys.argv)                                               
if(narg == 4):
    xbox = int(sys.argv[1])
    ybox = int(sys.argv[2])
    zbox = int(sys.argv[3])
else:                                                              
    sys.exit('3 arguments to be passed: ix, iy, iz')
############################# 

n_targ=0.000246356087788
b_targ = 1.22

print "## Target ##"
print "b=",b_targ, " n=",n_targ
print "  "

#Volume which was used to computed N(Mh)  (=n(Mh) * V)  #VGP:???
vol=500**3

#------------------------
fix_fsat = True   ; fsat0 = 0.2
fix_mu =   False  ; mu0 = 12.
doplot =   False
#-----------------------

# Fixed parameters 
sig = 0.12
alpha = 0.8

# OuterRim file
halodir = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/'
outdir = halodir+'ascii/OuterRim_STEP266_z0.865/subvols27/'
root = outdir+'OuterRim_STEP'+str(istep)+'_fofproperties_'
file = root+str(xbox)+str(ybox)+str(zbox)+'.txt' 

# Read Mh
data = np.genfromtxt(file)
logMh = data[:,6] #Halo masses
print logMh ; sys.exit()

## Read HMF ------ to be updated to reading the new tables when ready 
#Mh, N(Mh), b(Mh)
File = '/users/savila/ELGs_eBOSS/CUTE_box/plots/Halos/Mh_n_b_Halostest_halfres.txt'
array = np.loadtxt(File)
hmf = array[:, 1]
bh = array[:, 2]
mhist = array[:, 0]
#-------------------------------------------------------

#Grid for mu and fsat
if fix_fsat:
    fsat_arr = np.array([fsat0])  
else:
    fsat_arr = np.array([20.,0.2]) #np.arange(0.002, 0.3 ,0.005) 

if fix_mu:
    mu_arr = np.array([mu0])
else:
    mu_arr = np.arange(11.0,15.0,0.01)
print type(fsat_arr),type(mu_arr)

# Set grid arrays
bgal= np.zeros((len(fsat_arr),len(mu_arr)))
delt = np.zeros_like(bgal)
Ac = np.zeros_like(bgal)
As = np.zeros_like(bgal)

#Loop through the paramters (fsat and mu)
#and compute (n,b) 
for ii,fsat in enumerate(fsat_arr):
    for jj,mu in enumerate(mu_arr):
        logM0 = mu - 0.1
        logM1 = mu + 0.3

        # Calculate number of centrals
        Integrand = hmf*gcent(mhist, mu, sig)
        Ic = integrate.simps(Integrand)
        ncen = n_targ*(1.-fsat)*Ic
        Ac[ii,jj] = ncen/Ic
        if (ncen<=0.):
            print ('ERROR no centrals for fsat,mu=',fsat,mu)
            bgal[ii,jj] = 999. ; delt[ii,jj] =  999.
            break

        # Calculate number of satellites
        As[ii,jj] = 0. ; nsat = 0.
        if (fsat>0.):
            Integrand = hmf*gsat(mhist, logM0, logM1, alpha)
            Is = integrate.simps(Integrand)
            nsat = n_targ*fsat*Is
            As[ii,jj] = nsat/Is

        # Calculate the galaxy bias
        Integrand = bh*hmf*gcent(mhist, mu, sig)
        Ic = integrate.simps(Integrand)
        bcen = Ac[ii,jj]*Ic

        bsat = 0.
        if (fsat>0.):
            Integrand = bh*hmf*gsat(mhist, logM0, logM1, alpha)
            Is = integrate.simps(Integrand)
            bsat = As[ii,jj]*Is

        bgal[ii,jj] = (bcen+bsat)/(nsat+ncen)

        #Compute the "error" wrt the target quantities
        delt[ii,jj] = delta_sq(bgal[ii,jj], b_targ)

#Find the parameters that give us the minimum "error"
print "delt=",delt
m,n= np.where(delt==np.min(delt))
print('Mimimum:')
k=m
l=n
print("Input b, n=",b_targ, n_targ
print("Delta=",delt[k,l][0]," b=",bgal[k,l][0],\
          "fsat=",fsat_arr[k][0],"mu=",mu_arr[l][0])
print("Fix params: alpha= ",alpha," sigma= ",sig)
print("Derived params: Ac= ",Ac[k,l][0]," As= ",As[k,l][0])
print("Re-scaled params: logM0=",mu_arr[l][0] - 0.1," logM1=",mu_arr[l][0] + 0.3)

if doplot:
    plotname = outdir+'/plots/bestfit_gauss.pdf'
    fig = plt.figure(figsize = (8., 9.))
    xtit = '${\\rm logM_{h}}$' ; ytit = '${\\rm logN}$'
    plt.xlabel(xtit) ; plt.ylabel(ytit)
    plt.ylim([-3, 3])

    Integrand = hmf*gcent(mhist, mu_arr[l][0], sig)
    Ic = integrate.simps(Integrand)
    ncen = n_targ*(1-fsat_arr[k][0])*Ic
    plt.plot(mhist, np.log10(ncen), label='${\\rm N_{cen}}$') 
    
    if (fsat_arr[k][0]>0.):
        Integrand = hmf*gsat(mhist, mu_arr[l][0] - 0.1, mu_arr[l][0] + 0.3, alpha)
        Is = integrate.simps(Integrand)
        nsat = n_targ*fsat_arr[k][0]*Is
        plt.plot(mhist, np.log10(nsat), label='${\\rm N_{sat}}$')
        
    leg = plt.legend(loc=1)
    leg.draw_frame(False)
    plt.show()
    fig.savefig(plotname)
    print ('Ouput: ', plotname)
