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

# Create the files including a header
icount = 0
for ix in range(ncell):
    for iy in range(ncell):
        for iz in range(ncell):
            icount += 1
            filename = root+str(ix)+str(iy)+str(iz)+'.txt'
            # Create the file with a header
            with open(filename, 'w') as outf:
                outf.write('# fof_halo_center_x (Mpc/h), fof_halo_center_y (Mpc/h), fof_halo_center_z (Mpc/h), log10(fof_halo_mass /Msun/h), fof_halo_tag \n')
                outf.closed

## Read each galaxy and write it in the correct file
nroot = halodir+'HaloCatalog/STEP'+str(istep)  
files = glob.glob(nroot+'/*'+str(istep)+'*.fofproperties#*')
for infile in files:
    # Read the file
    tag = gio.gio_read(infile,"fof_halo_tag")
    mass = gio.gio_read(infile,"fof_halo_mass")
    xc = gio.gio_read(infile,"fof_halo_center_x")
    yc = gio.gio_read(infile,"fof_halo_center_y")
    zc = gio.gio_read(infile,"fof_halo_center_z")
    
    size = len(xc)
    for i,x in enumerate(xc):
        y = yc[i] ; z = zc[i]
        ix = int(x*float(ncell)/lbox) #; print ix
        iy = int(y*float(ncell)/lbox) #; print iy
        iz = int(z*float(ncell)/lbox) #; print iz
        ioc = ix + ncell*iy +ncell2*iz
        outfile = root+str(ix)+str(iy)+str(iz)+'.txt'

        if (ioc>ncell3 or ioc<0):
            print('ERROR: Cell ',ioc,' of ',ncell3)
            print('ERROR: Cell ',ix,iy,iz)
            print('ERROR: Position ',x,y,z)
            sys.exit()
        if(not os.path.isfile(outfile)):
            print('ERROR: Not found ',outfile)
            sys.exit()

        if (mass[i]>0.):
            lmass = np.log10(mass[i])
            #tofile = (x,y,z,lmass,tag[i])
            tofile = (x,y,z,lmass)

            with open(outfile, 'a') as outf:
#                np.savetxt(outf, tofile, fmt=('%10.5f %10.5f %10.5f %10.5f %i'))
                np.savetxt(outf, tofile, fmt=('%10.5f %10.5f %10.5f %10.5f'))
                outf.closed
        sys.exit()
print ('Output in ',root)
