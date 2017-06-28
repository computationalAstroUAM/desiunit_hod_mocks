import sys, os, glob, getpass
sys.path.append('/mnt/lustre/eboss/OuterRim/genericio/python/')
import genericio as gio
import numpy as np
import random
import matplotlib
matplotlib.use('Agg') 
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#############################                                      
#                                                                  
#  Input ARGUMENTS                                                 
#                                                                  
narg = len(sys.argv)                                               
if(narg == 3):
    istep = int(sys.argv[1])                                     
    lcut = sys.argv[2]
else:                                                              
    sys.exit('2 arguments to be passed: step, lcut')

############################# 

# Directory with the OuterRim simulation haloes
halodir = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/'

# OuterRim simulation characteristics (FOF b=0.168 here)
mp  = 1.9E+09 # Msol/h
lbox= 3000.   # Mpc/h

# Get the conversion between the name of the time step and redshift
step = np.genfromtxt(halodir+'step_redshift.txt',usecols=0,dtype=int)
redshift = np.genfromtxt(halodir+'step_redshift.txt',usecols=1)     

zz = redshift[np.where(step == istep)]
outdir = halodir+'ascii/vol'+lcut+'Mpc/OuterRim_STEP'+\
    str(istep)+'_z'+str(zz[0])+'/'    

if not os.path.exists(outdir):
    os.makedirs(outdir)
print 'Processing redshift ',zz

######## FOF properties
fofproperties = True
if (fofproperties):
    print 'Halo catalogue --------------------'

    #Output file 
    outstep = 500000
    outfile = outdir+'OuterRim_STEP'+str(istep)+'_fofproperties_'+lcut+'Mpc.txt'
    if os.path.isfile(outfile):
        os.remove(outfile)

    # File
    nroot = halodir+'HaloCatalog/STEP'+str(istep)    
    files = glob.glob(nroot+'/*'+str(istep)+'*.fofproperties#*')
    first = True
    for infile in files:
        # Read the file
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
        size = len(vz)

        ind = np.where((xc>0.) & (xc<float(lcut)) & \
                           (yc>0.) & (yc<float(lcut)) & \
                           (zc>0.) & (zc<float(lcut)) )
        if (np.shape(ind)[1] > 0):
            print np.shape(ind)[1],' haloes in ',infile
        
            tofile = zip(count[ind],tag[ind],mass[ind],\
                             xc[ind],yc[ind],zc[ind],\
                             xm[ind],ym[ind],zm[ind],\
                         vx[ind],vy[ind],vz[ind])

            if (first):
                first = False
                print '----------------------------'
                gio.gio_inspect(infile)
                print '----------------------------'
            
                # Write header
                with open(outfile, 'w') as outf:  
                    outf.write('# fof_halo_count, fof_halo_tag, fof_halo_mass, fof_halo_center_x, fof_halo_center_y, fof_halo_center_z, fof_halo_mean_x, fof_halo_mean_y, fof_halo_mean_z, fof_halo_mean_vx, fof_halo_mean_vy, fof_halo_mean_vz \n')
                    np.savetxt(outf,tofile,fmt=('%i %i %3.5e %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f'))    
                    outf.closed 
            else:
                with open(outfile, 'a') as outf:                            
                    np.savetxt(outf,tofile,fmt=('%i %i %3.5e %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f'))    
                    outf.closed 
            #print outfile ; sys.exit()
    print '----------------------------'
    print 'Output fof properties: ',outfile
                
    # Test the output file
    x,y,z = np.loadtxt(outfile,usecols=(3,4,5),unpack=True)
    print 'x (min,max), y (min,max), z(min,max):'
    print min(x),max(x),min(y),max(y),min(z),max(z)

    # Downsample to plot
    num = len(x) ; print num
    idx = np.arange(num) ; random.shuffle(idx)
    val = num/100

    # Plot
    fig = plt.figure() ; ax = fig.add_subplot(111,projection='3d')
    xtit ='x (Mpc/h)' ; ytit ='y (Mpc/h)'; ztit ='z (Mpc/h)'
    ax.set_xlabel(xtit) ; ax.set_ylabel(ytit) ; ax.set_zlabel(ztit)
    ax.scatter(xs=x[idx[:val]],ys=y[idx[:val]],zs=z[idx[:val]])

    plotfile = outdir+'OuterRim_STEP'+str(istep)+'_xyz_'+lcut+'Mpc.pdf'
    fig.savefig(plotfile)
    print 'Test plot: ',plotfile

######## Particles
particles = True
if (particles):
    print 'Particle catalogue --------------------'

    #Output file 
    outfile = outdir+'OuterRim_STEP'+str(istep)+'_particles_'+lcut+'Mpc.txt'
    if os.path.isfile(outfile):
        os.remove(outfile)

    # File
    nroot = halodir+'HaloCatalog/STEP'+str(istep)    
    files = glob.glob(nroot+'/*'+str(istep)+'*.haloparticles#*')
    first = True
    for infile in files:
        # Read the file
        xc = gio.gio_read(infile,"x")
        yc = gio.gio_read(infile,"y")
        zc = gio.gio_read(infile,"z")
        vx = gio.gio_read(infile,"vx")
        vy = gio.gio_read(infile,"vy")
        vz = gio.gio_read(infile,"vz")
        ids= gio.gio_read(infile,"id")
        tag= gio.gio_read(infile,"fof_halo_tag")

        # Output
        size = len(xc)

        ind = np.where((xc>0.) & (xc<float(lcut)) & \
                           (yc>0.) & (yc<float(lcut)) & \
                           (zc>0.) & (zc<float(lcut)) )
        if (np.shape(ind)[1] > 0):
            print np.shape(ind),' haloes in ',infile
            
            tofile = zip(xc[ind],yc[ind],zc[ind],\
                             vx[ind],vy[ind],vz[ind],\
                             ids[ind],tag[ind])
            if (first):
                first = False
                gio.gio_inspect(infile)
            
                # Write header
                with open(outfile, 'w') as outf:
                    outf.write('# x (cMpc/h), y (cMpc/h), z (cMpc/h), vx (ckm/s), vy (ckm/s), vz (ckm/s), id (from 0 to the total number of particles), fof_halo_tag (to be matched to the halo file) \n') 
                    np.savetxt(outf,tofile,fmt=('%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %i %i'))    
                    outf.closed 
            else:
                with open(outfile, 'a') as outf:                            
                    np.savetxt(outf,tofile,fmt=('%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %i %i'))   
                    outf.closed                     
            #print outfile ; sys.exit()
    print '----------------------------'
    print 'Output particles: ',outfile

