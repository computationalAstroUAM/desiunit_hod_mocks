import scipy.integrate as integrate
import numpy as np
import sys
from decimal import *
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

getcontext().prec = 6 #For 'decimal' precision

def Calcgsat(M, M0, M1, alpha):
    Nsat_Mh = 0. 
    if (M>M0):
        Nsat_Mh = (((10**M) - 10**M0)/10**M1)**alpha
    return Nsat_Mh
gsat = np.vectorize(Calcgsat)

def Calcgcent(M, mu, sig):
    return np.exp(-(M - mu)**2/(2.0*sig**2))/(sig*np.sqrt(np.pi*2.0))
gcent = np.vectorize(Calcgcent)

def CalcAsat(nd, vol, fsat, Is):
    return nd*vol*fsat/Is
Asat = np.vectorize(CalcAsat)

def CalcAcen(nd, vol, fsat, Ic):
    return nd*vol*(1.0 - fsat)/Ic
Acen = np.vectorize(CalcAcen)

def mu2m0(mu):
    return mu - 0.1
def mu2m1(mu):
    return mu + 0.3

#Metrics that we want to minimisize (or "error")
def delta_sq(b_mod, b_targ):
    delt =((b_mod**2 - b_targ**2)/b_targ**2 )**2
    return delt

#############################                                      
#                                                                  
#  Input ARGUMENTS                                                 
#                                                                  
n_targ=0.000246356087788
b_targ = 1.22

# Fixed parameters 
sig = 0.12
alpha = 1. #0.8
############################# 
fix_fsat = False  ; fsat0 = 0.2
fix_mu =   False  ; mu0 = 12.
doplot =   False
verbose =  False
############################# 

if verbose:
    print "## Target ##"
    print "b=",b_targ, " n=",n_targ
    print "  "

# Read HMF (number of haloes per dlogMh per volume)
halodir = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/'
mfile = halodir+'hmf.txt'
hmf1, mhmed1 = \
    np.loadtxt(mfile, usecols= (2,3), unpack=True )
mhist = mhmed1 #mhmed1[1:-1]
hmf = hmf1 #hmf1[1:-1]

# Read bias function
bfile = halodir+'bias_rl20.0_rh80.0.txt'
bh1 = np.loadtxt(bfile, usecols= (1,), unpack=True )
bh = bh1 #bh1[1:-1]

#-------------------------------------------------------

#Grid for mu and fsat
if fix_fsat:
    fsat_arr = np.array([fsat0])
else:
    fsat_arr = np.arange(0.002, 0.3 ,0.005)

if fix_mu:
    mu_arr = np.array([mu0])
else:
    mu_arr = np.arange(10.5,13.5,0.005)

if verbose:
    print ('Halo masses:',mhist)
    print ('fsat values:',fsat_arr)
    print ('mu values:',mu_arr)

# Set grid arrays
bgal= np.zeros((len(fsat_arr),len(mu_arr)))
delt = np.zeros_like(bgal)
Ac = np.zeros_like(bgal)
As = np.zeros_like(bgal)

#Loop through the paramters (fsat and mu)
#and compute (n,b) 
for ii,fsat in enumerate(fsat_arr):
    for jj,mu in enumerate(mu_arr):
        logM0 = mu2m0(mu)
        logM1 = mu2m1(mu)

        # Calculate Ac
        ncen = n_targ*(1.-fsat)
        if (ncen<=0.):
            print ('ERROR no centrals for fsat,mu=',fsat,mu)
            bgal[ii,jj] = 999. ; delt[ii,jj] =  999.
            break

        Integrand = hmf*gcent(mhist, mu, sig)
        Ic = integrate.simps(Integrand,mhist)
        Ac[ii,jj] = ncen/Ic #; print Ic,'# Ic' 

        # Calculate As
        nsat = n_targ*fsat 

        As[ii,jj] = 0. 
        if (fsat>0.):
            Integrand = hmf*gsat(mhist, logM0, logM1, alpha)
            Is = integrate.simps(Integrand,mhist)
            if (Is>0.):
                As[ii,jj] = nsat/Is

        # Calculate the galaxy bias
        Integrand = bh*hmf*gcent(mhist, mu, sig)
        Ic = integrate.simps(Integrand,mhist)
        Is = 0.
        if (fsat>0.):
            Integrand = bh*hmf*gsat(mhist, logM0, logM1, alpha)
            Is = integrate.simps(Integrand,mhist)

        bgal[ii,jj] = (Ac[ii,jj]*Ic + As[ii,jj]*Is)/n_targ

        #Compute the "error" wrt the target quantities
        delt[ii,jj] = delta_sq(bgal[ii,jj], b_targ)

#Find the parameters that give us the minimum "error"
if verbose:
    print("  ")
    print("delt=",delt) 
    print("  ")
print('Input b= {}, n= {}'.format(b_targ, n_targ))
print('Fix params: alpha= {}, sigma= {}'.format(alpha,sig))
m,n= np.where(delt==np.nanmin(delt))
k=m
l=n 
print('Minimum delta= {}'.format(delt[k,l][0]))
print('  ')
print('-----Best fit-----')
b_bf = bgal[k,l][0] ; fsat_bf = fsat_arr[k][0] ; mu_bf = mu_arr[l][0]
print('b= {:10.4f}, fsat= {:10.4f}, mu= {:10.4f}'.format(b_bf,fsat_bf,mu_bf))
ac_bf = Ac[k,l][0] ; as_bf = As[k,l][0] 
print('Derived params: Ac= {}, As= {}'.format(ac_bf,as_bf))
m0_bf = mu2m0(mu_bf) ; m1_bf = mu2m1(mu_bf)
print('Re-scaled params: logM0= {:10.4f}, logM1= {:10.4f}'.format(m0_bf,m1_bf))
print('  ')
print('-----Test-----')
Integrand = hmf*gcent(mhist, mu_bf, sig)
Ic = integrate.simps(Integrand,mhist)
Integrand = hmf*gsat(mhist, m0_bf, m1_bf, alpha)
Is = integrate.simps(Integrand,mhist)
ncs = ac_bf*Ic + as_bf*Is
print('ngal(As,Ac)= {}, ngal/n_targ= {}'.format(ncs,ncs/n_targ))

Integrand = bh*hmf*gcent(mhist, mu_bf, sig)
Ic = integrate.simps(Integrand,mhist)
Integrand = bh*hmf*gsat(mhist, m0_bf, m1_bf, alpha)
Is = integrate.simps(Integrand,mhist)
bcs = (ac_bf*Ic + as_bf*Is)/n_targ
print('bias(As,Ac)= {:10.4f}, bias/b_targ= {:10.4f}'.format(bcs,(bcs/b_targ)))
print('--------------')

if doplot:
    plotname = halodir+'plots/bestfit_gauss.png'
    fig = plt.figure(figsize = (8., 9.))
    xtit = '${\\rm logM_{h}}$' ; ytit = '${\\rm logN}$'
    plt.xlabel(xtit) ; plt.ylabel(ytit)
    plt.ylim([-3, 3])

    ncen_mh = Ac[k,l][0]*gcent(mhist, mu_arr[l][0], sig)
    plt.plot(mhist, np.log10(ncen_mh), label='${\\rm N_{cen}}$') 
    
    if (fsat_arr[k][0]>0.):
        nsat_mh = As[k,l][0]*gsat(mhist, mu_arr[l][0] - 0.1, mu_arr[l][0] + 0.3, alpha)
        plt.plot(mhist, np.log10(nsat_mh), label='${\\rm N_{sat}}$')
        
    leg = plt.legend(loc=1)
    leg.draw_frame(False)
    plt.show()
    fig.savefig(plotname)
    print ('Ouput: ', plotname)
