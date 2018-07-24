import sys, os, glob, getpass
#sys.path.append('/mnt/lustre/eboss/genericio/python/')
#sys.path.append('../tools')
import genericio as gio
import numpy as np
import matplotlib ; matplotlib.use('Agg')                                    
from matplotlib import pyplot as plt  
from distinct_colours import get_distinct 
cols = get_distinct(10) 

# Directory with the OuterRim simulation haloes
halodir = '/mnt/lustre/eboss/OuterRim/'

# OuterRim simulation characteristics (FOF b=0.168 here)
mp  = 1.9E+09 # Msol/h
lbox= 3000.   # Mpc/h

# Get the conversion between the name of the time step and redshift
step = np.genfromtxt(halodir+'step_redshift.txt',usecols=0,dtype=int)
redshift = np.genfromtxt(halodir+'step_redshift.txt',usecols=1)

# Initialize the parameters for the figure ------------------
plt.rcParams['legend.numpoints'] = 1 
plt.rcParams['axes.labelsize'] = 10.0 ; fs = 20 
plt.rcParams['lines.linewidth'] = 2 
fig = plt.figure(figsize=(8.5,9.))

xtit = "${\\rm log}_{10}(\\rm{M/M_{\odot}}h^{-1})$"
ytit = "${\\rm log}_{10}(\Phi/ Mpc^{-3}h^3 {\\rm dlog}_{10}M)$"
                                                                             
xmin = 10. ; xmax = 16.                                                      
ymin = -6.5 ; ymax = 0.    

jj = 111                                                                     
ax = fig.add_subplot(jj) 
ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax)                              
ax.set_xlabel(xtit,fontsize = fs) ; ax.set_ylabel(ytit,fontsize = fs)        

# Generate halo mass functions as a test
mmin = 10. ; mmax = 16. ; dm = 0.1
mbins = np.arange(mmin,mmax,dm)                                              
mhist = mbins + dm*0.5  
#-------------------------------------------------------------                

#for iz,istep in enumerate(step):  # Loop over all redshifts
# Loop over a subset of redshifts
for iz,istep in enumerate([499, 307, 279, 247, 224]):    
    zz = redshift[np.where(step == istep)]
    print 'Processing snapshot at redshift ',zz
    nroot = halodir+'HaloCatalog/STEP'+str(istep)    

    # Initialize to 0 the halo mass functions
    ycount, yh = [np.zeros(len(mhist)) for _ in range(2)]

    # Loop over each of the sub volumes the 
    for inf,infile in enumerate(glob.glob(nroot+'/*'+str(istep)+'*#*')) :
        # Print out once the information stored in each file
        if (iz == 0 and inf == 0):
            gio.gio_inspect(infile)
            print '----------------------------'
 
        # Read the number of particles per halo
        in_count = mp*gio.gio_read(infile, "fof_halo_count")    
        ind = np.where(in_count >0.) ; count = np.log10(in_count[ind])    
        H, bins_edges = np.histogram(count,bins=np.append(mbins,mmax))
        ycount = ycount + H

        # FOF mass (Msun/h)
        in_mh = gio.gio_read(infile, "fof_halo_mass")
        ind = np.where(in_mh >0.) ; mh = np.log10(in_mh[ind])
        H, bins_edges = np.histogram(count,bins=np.append(mbins,mmax))
        yh = yh + H
        
    ycount = ycount/dm/(lbox**3)
    yh = yh/dm/(lbox**3)

    ind = np.where(ycount >0.)
    ax.plot(mhist[ind],np.log10(ycount[ind]),\
                color=cols[iz],label='z='+str(zz))

    ind = np.where(yh >0.)
    ax.plot(mhist[ind],np.log10(yh[ind]),\
                color=cols[len(cols)-iz-1],linestyle=':',label=' from mass')


# Legend
plt.legend(loc=3,prop={'size':(fs-2)}) 


# Directory with outputs (it'll be generated if it doesn't exist)
outdir = '/users/'+getpass.getuser()+'/Outputs/out_shams/'
if (not os.path.exists(outdir)): os.makedirs(outdir)

# Save figure
plotfile = outdir+'outerrim_mf.pdf'
fig.savefig(plotfile)
print 'Output: ',plotfile
