import sys, os, glob, getpass
#sys.path.append('/mnt/lustre/eboss/genericio/python/')
#sys.path.append('../tools')
import genericio as gio
import numpy as np

# Directory with the OuterRim simulation haloes
halodir = '/mnt/lustre/eboss/OuterRim/'

# OuterRim simulation characteristics (FOF b=0.168 here)
mp  = 1.9E+09 # Msol/h
lbox= 3000.   # Mpc/h

# Get the conversion between the name of the time step and redshift
step = np.genfromtxt(halodir+'step_redshift.txt',usecols=0,dtype=int)
redshift = np.genfromtxt(halodir+'step_redshift.txt',usecols=1)     

# Loop over a subset of redshifts
for iz,istep in enumerate([203, 266, 300]):    
    zz = redshift[np.where(step == istep)]
    print 'Processing snapshot at redshift ',zz
    nroot = halodir+'HaloCatalog/STEP'+str(istep)    

    outdir = halodir+'ascii/OuterRim_STEP'+str(istep)+'_z'+str(zz[0])+'/'    
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Loop over each of the sub volumes the 
    for inf,infile in enumerate(glob.glob(nroot+'/*'+str(istep)+'*#*')) :
        # Print out (only once) the information stored in each file
        if (iz == 0 and inf == 0):
            gio.gio_inspect(infile)
            print '----------------------------'

        # Read the files
        count = gio.gio_read(infile,"fof_halo_count")
        tag = gio.gio_read(infile,"fof_halo_tag")
        mass = gio.gio_read(infile,"fof_halo_mass")
        xc = gio.gio_read(infile,"fof_halo_center_x")
        yc = gio.gio_read(infile,"fof_halo_center_y")
        zc = gio.gio_read(infile,"fof_halo_center_z")
        xm = gio.gio_read(infile,"fof_halo_mean_x")
        ym = gio.gio_read(infile,"fof_halo_mean_y")
        zm = gio.gio_read(infile,"fof_halo_mean_z")
        vx = gio.gio_read(infile,"fof_halo_mean_vx")
        vy = gio.gio_read(infile,"fof_halo_mean_vy")
        vz = gio.gio_read(infile,"fof_halo_mean_vz")

        # Output
        tofile = zip(count,tag,mass,xc,yc,zc,xm,ym,zm,vx,vy,vz)

        outfile = outdir+'fofproperties_'+str(inf)+'.txt'
        with open(outfile, 'w') as outf:                            
            np.savetxt(outf,tofile,fmt=('%i %i %3.5e %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f'))    
            outf.closed 
print 'Output: ',outfile
