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

#def Calcgcent(M, mu, sig):
#    return np.exp(-(M - mu)**2/(2.0*sig**2))/(sig*np.sqrt(np.pi*2.0))

def CalcNcent(M, mu, sig,gamma):
	if (M<mu):
		return 1.0/(sig*np.sqrt(np.pi*2.0))*np.exp(-(M-mu)**2/(2.0*sig**2))
	else:
		return 1.0/(sig*np.sqrt(np.pi*2.0))*(10**((M-mu)*gamma))
gcent = np.vectorize(CalcNcent)


def CalcAsat(nd, vol, fsat, Is):
    return nd*vol*fsat/Is
Asat = np.vectorize(CalcAsat)

def CalcAcen(nd, vol, fsat, Ic):
    return nd*vol*(1.0 - fsat)/Ic
Acen = np.vectorize(CalcAcen)

def mu2m0(mu):
    return mu - 0.05
def mu2m1(mu):
    return mu + 0.35

#Metrics that we want to minimisize (or "error")
def delta_sq(b_mod, b_targ):
    delt =((b_mod**2 - b_targ**2)/b_targ**2 )**2
    return delt

#############################                                      
#                                                                  
#  Input ARGUMENTS                                                 
#    

target=1
corr=1
times7=True
if (target == 0):
	#V3 
	print("\n** Target V3 **\n")
	#n_targ = 0.000246549484522
	#b_targ = 1.3193
elif (target == 1):
	#V7
	print("\n** Target V7 **\n")
	if corr==0:
		n_targ = 0.00021873633933254137 ## mean of SGC+NGC
	else:
		n_targ = 0.00021873633933254137**2/0.00021292622222222226 # corrected (integral biased)
	b_targ = 1.3198 ## corrected for Kaiser
	


# Fixed parameters 
sig = 0.08
alpha = 0.9
gamma = -1.4
############################# 
fix_fsat = False  ; fsat0 = 0.15
fix_mu =   False  ; mu0 = 12.
doplot =   False
verbose =  True
############################# 

if verbose:
    print("## Target ##")
    print("b=",b_targ, " n=",n_targ)
    print("  ")

# Read HMF (number of haloes per dlogMh per volume)
#halodir = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/'
#mfile = halodir+'hmf.txt'
mfile = '/users/savila/ELGs_eBOSS/python_scripts/hmf_finest.txt'
hmf1, mhmed1 = \
    np.loadtxt(mfile, usecols= (2,3), unpack=True )
mhist = mhmed1 #mhmed1[1:-1]
hmf = hmf1 #hmf1[1:-1]

# Read bias function
bfile = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/bias_rl20.0_rh80.0.txt'
M,b = np.loadtxt(bfile, usecols= (0,1), unpack=True )
p4 = np.poly1d(np.polyfit(M[:76], b[:76], 4))

bh = p4(mhmed1)

#-------------------------------------------------------

#Grid for mu and fsat
if fix_fsat:
    fsat_arr = np.array([fsat0])
else:
    fsat_arr = np.arange(0.00, 0.81,0.01)

if fix_mu:
    mu_arr = np.array([mu0])
else:
    mu_arr = np.arange(10.5,13.5,0.001)

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
#out = np.array([])
out = []
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

        Integrand = hmf*gcent(mhist, mu, sig,gamma)
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
        Integrand = bh*hmf*gcent(mhist, mu, sig,gamma)
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
    n= np.where(delt[ii,:]==np.nanmin(delt[ii,:]))
    l=n[0][0] 
    k=ii
#    print('Minimum delta= {}'.format(delt[ii,l]))
    print('  ')
    print('-----Best fit-----')
    b_bf = bgal[ii,l] 
    fsat_bf = fsat_arr[ii] 
    mu_bf = mu_arr[l]
    if verbose:
        print('b= {:10.4f}, fsat= {:10.4f}, mu= {:10.4f}'.format(b_bf,fsat_bf,mu_bf))
    ac_bf = Ac[k,l] ; as_bf = As[k,l] 
    if verbose:
        print('Derived params: Ac= {}, As= {}'.format(ac_bf,as_bf))
    m0_bf = mu2m0(mu_bf) ; m1_bf = mu2m1(mu_bf)
#    out = np.append(out,np.array([[fsat,mu_bf,ac_bf,as_bf]])) 
    if (times7):
    	out.append([fsat,mu_bf,ac_bf*7,as_bf*7])
    else:
    	out.append([fsat,mu_bf,ac_bf,as_bf])
    if verbose:
        print('Re-scaled params: logM0= {:10.4f}, logM1= {:10.4f}'.format(m0_bf,m1_bf))
        print('  ')
    if verbose:
        print('-----Test-----')
    Integrand = hmf*gcent(mhist, mu_bf, sig,gamma)
    Ic = integrate.simps(Integrand,mhist)
    Integrand = hmf*gsat(mhist, m0_bf, m1_bf, alpha)
    Is = integrate.simps(Integrand,mhist)
    ncs = ac_bf*Ic + as_bf*Is
    if verbose:
        print('ngal(As,Ac)= {}, ngal/n_targ= {}'.format(ncs,ncs/n_targ))

    Integrand = bh*hmf*gcent(mhist, mu_bf, sig,gamma)
    Ic = integrate.simps(Integrand,mhist)
    Integrand = bh*hmf*gsat(mhist, m0_bf, m1_bf, alpha)
    Is = integrate.simps(Integrand,mhist)
    bcs = (ac_bf*Ic + as_bf*Is)/n_targ
    if verbose:
        print('bias(As,Ac)= {:10.4f}, bias/b_targ= {:10.4f}'.format(bcs,(bcs/b_targ)))
    print('For fsat={:1.2f} Use mu,Ac,As={:1.8f} {:1.10f} {:1.10f}'.format(fsat,mu_bf,ac_bf,as_bf))
    print('--------------')
print("out",out)
if (times7):
	if (corr==0):
		np.savetxt('HOD_gaussPL_PL_x7.params',out,fmt='%.2f %.6f %e %e')
	else:
		np.savetxt('HOD_gaussPL_PL_corr_x7.params',out,fmt='%.2f %.6f %e %e')
else:
	if (corr==0):
		np.savetxt('HOD_gaussPL_PL.params',out,fmt='%.2f %.6f %e %e')
	else:
		np.savetxt('HOD_gaussPL_PL_corr.params',out,fmt='%.2f %.6f %e %e')
	
