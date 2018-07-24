import sys, os, glob
#sys.path.append('/mnt/lustre/eboss/genericio/python/')
import genericio as gio
import numpy as np

# Directory with the OuterRim simulation haloes
halodir = '/mnt/lustre/eboss/OuterRim/'

# Get the conversion between the name of the time step and redshift
step = np.genfromtxt(halodir+'step_redshift.txt',usecols=0,dtype=int)
redshift = np.genfromtxt(halodir+'step_redshift.txt',usecols=1)

istep = step[2]
nroot = halodir+'HaloCatalog/STEP'+str(istep)    

files = glob.glob(nroot+'/*'+str(istep)+'*#*')
infile = files[0]

gio.gio_inspect(infile)

x = gio.gio_read(infile, "fof_halo_mean_x")    
print x
