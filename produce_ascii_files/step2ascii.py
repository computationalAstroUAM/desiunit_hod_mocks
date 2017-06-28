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
outfile = outdir+'OuterRim_STEP'+str(istep)+'_fofproperties_'+str(ifile)+'.txt'
if os.path.isfile(outfile):
    os.remove(outfile)

# File
nroot = halodir+'HaloCatalog/STEP'+str(istep)    
files = glob.glob(nroot+'/*'+str(istep)+'*.fofproperties#'+ifile)
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

    for i in range(0,size,outstep):
        j = min(i+outstep-1,size-1)
        #print i, j, size-1
        tofile = zip(count[i:j],tag[i:j],mass[i:j],\
                         xc[i:j],yc[i:j],zc[i:j],\
                         xm[i:j],ym[i:j],zm[i:j],\
                         vx[i:j],vy[i:j],vz[i:j])

        if (i == 0):
            with open(outfile, 'w') as outf:   
                outf.write('# fof_halo_count, fof_halo_tag, fof_halo_mass, fof_halo_center_x, fof_halo_center_y, fof_halo_center_z, fof_halo_mean_x, fof_halo_mean_y, fof_halo_mean_z, fof_halo_mean_vx, fof_halo_mean_vy, fof_halo_mean_vz \n')
                np.savetxt(outf,tofile,fmt=('%i %i %3.5e %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f'))    
                outf.closed 
        else:
            with open(outfile, 'a') as outf:                            
                np.savetxt(outf,tofile,fmt=('%i %i %3.5e %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f'))    
                outf.closed 
                
