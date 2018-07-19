import scipy.integrate as integrate
import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

def Calcgsat(M, M0, M1, alpha):
        Nsat_Mh = (((10**M) - 10**M0)/10**M1)**alpha
        if (M>M0):
                return Nsat_Mh
        else:
                return 0.

gsat = np.vectorize(Calcgsat)

def Calcgcent(M, mu, sig):
        return 1.0/(sig*np.sqrt(np.pi*2.0))*np.exp(-(M - mu)**2/(2.0*sig**2))

gcent = np.vectorize(Calcgcent)

def CalcAsat(n, V, fsat, Is):
        return n*V*fsat/Is

Asat = np.vectorize(CalcAsat)

def CalcAcen(n, V, fsat, Ic):
        return n*V*(1.0 - fsat)/Ic

Acen = np.vectorize(CalcAcen)

doplot = True

#Input variables

data = np.genfromtxt("/mnt/lustre/eboss/OuterRim/OuterRim_sim/ascii/vol200Mpc/OuterRim_STEP266_z0.865/OuterRim_STEP266_fofproperties_200Mpc.txt")
M = data[:,2] #Halo masses
logM = np.log10(M)

V=200**3
n=0.000246356087788
fsat=0.2

sig = 0.12
mu = 12.0

logM0 = mu
alpha = 0.8
logM1 = 1.0 + mu

# Set up histogram values
logMmin = 8.0 #Minimum log(mass)
logMmax = 16.0 #Maximum log(mass)
dlogM = 0.1 #Histogram bin width

mbins = np.arange(logMmin,logMmax,dlogM) #Creates the bins
mhist = mbins + dlogM*0.5 #Defines the middle of each bin

H, bins_edges = np.histogram(logM, bins=np.append(mbins, logMmax)) #Creates a histogram, used to define dNh/dlogMh
dNhdlogMh = H

Integrand1 = dNhdlogMh*gsat(mhist, logM0, logM1, alpha)
Integrand2 = dNhdlogMh*gcent(mhist, mu, sig)
Is = integrate.simps(Integrand1)
Ic = integrate.simps(Integrand2)

print ('Asat=', Asat(n, V, fsat, Is))
print ('Acen=', Acen(n, V, fsat, Ic))

if doplot:
	path2plot = "/users/ghthomas/output_plots/"
	plotname = "logMh vs logNcen and logNsat_200Mpc redux.png"
	fig = plt.figure(figsize = (8., 9.))
	xtit = '${\\rm logM_{h}}$'
	ytit = '${\\rm logN}$'
	plt.xlabel(xtit)
	plt.ylabel(ytit)
	plt.ylim([-3, 3])
	plt.plot(mhist, np.log10(Acen(n, V, fsat, Ic)*gcent(mhist, mu, sig)), label='${\\rm N_{cen}}$') 
	plt.plot(mhist, np.log10(Asat(n, V, fsat, Is)*gsat(mhist, logM0, logM1, alpha)), label='${\\rm N_{sat}}$')
	leg = plt.legend(loc=1)
	leg.draw_frame(False)
	plt.show()
	fig.savefig(path2plot + plotname)
	print ('Ouput: ', path2plot + plotname)
