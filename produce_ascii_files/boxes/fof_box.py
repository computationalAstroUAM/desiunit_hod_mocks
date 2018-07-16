import sys, os, glob, getpass
sys.path.append('/mnt/lustre/eboss/OuterRim/genericio/python/')
import genericio as gio
import numpy as np

# Directory with the OuterRim simulation haloes
halodir = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/'
istep = 266

outdir = halodir+'ascii/OuterRim_STEP'+str(istep)+'_z0.865/subvols27/'

root = outdir+'OuterRim_STEP'+str(istep)+'_fofproperties_'

# Subvolumes
lbox = 3000.
cell = lbox/3. ; print('Cell size (Mpc/h) =',cell)
ncell = int(lbox/cell) ; ncell2 = ncell*ncell
ncell3 = ncell2*ncell #; print ncell, ncell2, ncell3

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
# Minimum values
allx0, ally0, allz0 = [np.zeros(ncell) for i in range(3)]
for i in range(3):
    allx0[i] = i*cell
    ally0[i] = i*cell
    allz0[i] = i*cell
          
x0 = allx0[xbox] ; xL = x0 + cell
y0 = ally0[ybox] ; yL = y0 + cell
z0 = allz0[zbox] ; zL = z0 + cell
print type(x0), type(xL)
## Read each galaxy and write it in the correct file
nroot = halodir+'HaloCatalog/STEP'+str(istep)  
files = glob.glob(nroot+'/*'+str(istep)+'*.fofproperties#*')
ifile = 0
for infile in files:
    # Read the file
    tag = gio.gio_read(infile,"fof_halo_tag")
    mass = gio.gio_read(infile,"fof_halo_mass")
    xc = gio.gio_read(infile,"fof_halo_center_x")
    yc = gio.gio_read(infile,"fof_halo_center_y")
    zc = gio.gio_read(infile,"fof_halo_center_z")
    vx = gio.gio_read(infile,"fof_halo_mean_vx")
    vy = gio.gio_read(infile,"fof_halo_mean_vy")
    vz = gio.gio_read(infile,"fof_halo_mean_vz")

    # Output
    ind = np.where((mass>0.) & \
                       (xc>x0) & (xc<xL) & \
                       (yc>y0) & (yc<yL) & \
                       (zc>z0) & (zc<zL) )

    if (np.shape(ind)[1] > 0):
        lmass = np.log10(mass[ind])
        print lmass ; sys.exit()
        # Shift values so all the boxes start at (0,0,0)
        xs = x[ind] - x0
        ys = y[ind] - y0
        zs = z[ind] - z0

        tofile = zip(xs,ys,zs,\
                         vx[ind],vy[ind],vz[ind],\
                         lmass,tag[ind])

        if os.path.exists(filename):
            with open(outfile, 'a') as outf:
                outf.write(tofile+' \n') #; print outfile
                outf.closed
        else:
            with open(outfile, 'w') as outf:
                outf.write(tofile+' \n') #; print outfile
                outf.closed

    # Testing -------------
    ifile += 1
    if ifile>2:
        break
    #-------------------------

print ('Output in ',root)
