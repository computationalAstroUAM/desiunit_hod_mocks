import sys, os, glob, getpass
sys.path.append('/mnt/lustre/eboss/OuterRim/genericio/python/')
import genericio as gio
import numpy as np

#############################                                      
#                                                                  
#  Input ARGUMENTS                                                 
#                                                                  
narg = len(sys.argv)                                               
if(narg == 3):                                                     
    istep = int(sys.argv[1])                                     
    ifile = sys.argv[2]
else:                                                              
    sys.exit('2 arguments to be passed')                           
                                                                   
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
outdir = halodir+'ascii/OuterRim_STEP'+str(istep)+'_z'+str(zz[0])+'/'    
if not os.path.exists(outdir):
    os.makedirs(outdir)
print 'Processing redshift ',zz

#Output file 
outstep = 500000
outfile = outdir+'OuterRim_STEP'+str(istep)+'_particles_'+str(ifile)+'.txt'
if os.path.isfile(outfile):
    os.remove(outfile)

# File
nroot = halodir+'HaloCatalog/STEP'+str(istep)    
files = glob.glob(nroot+'/*'+str(istep)+'*.haloparticles#'+ifile)
if (len(files) > 1):
    print files
    print 'STOP: More than one file for ifile ',ifile
    sys.exit()
else:
    infile = files[0]
    gio.gio_inspect(infile)
    print '----------------------------'
    print 'Output: ',outfile

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
    size = len(vz)

    for i in range(0,size,outstep):
        j = min(i+outstep-1,size-1)
        #print i, j, size-1
        tofile = zip(xc[i:j],yc[i:j],zc[i:j],\
                         vx[i:j],vy[i:j],vz[i:j],\
                         ids[i:j],tag[i:j])

        if (i == 0):
            with open(outfile, 'w') as outf:
                outf.write('# x (cMpc/h), y (cMpc/h), z (cMpc/h), vx (ckm/s), vy (ckm/s), vz (ckm/s), id (from 0 to the total number of particles), fof_halo_tag (to be matched to the halo file) \n') 
                np.savetxt(outf,tofile,fmt=('%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %i %i'))    
                outf.closed 
        else:
            with open(outfile, 'a') as outf:                            
                np.savetxt(outf,tofile,fmt=('%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %i %i'))    
                outf.closed 
                
