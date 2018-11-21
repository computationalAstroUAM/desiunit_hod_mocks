import sys, os
import numpy as np

istep = 300
iz = 0.656

#241 1.055  R
#253 0.959  R
#279 0.779  R 
#300 0.656  R  

# Directory with the OuterRim simulation haloes
halodir = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/'

# Create output directory if it does not exist
outdir = halodir+'ascii/OuterRim_STEP'+str(istep)+'_z'+str(iz)+'/subvols27/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

root = outdir+'OuterRim_STEP'+str(istep)+'_fofproperties_'

# Subvolumes
lbox = 3000.
cell = lbox/3. ; print('Cell size (Mpc/h) =',cell)
ncell = int(lbox/cell) ; ncell2 = ncell*ncell
ncell3 = ncell2*ncell #; print ncell, ncell2, ncell3

# Remove files if they already exist
for ix in range(ncell):
    for iy in range(ncell):
        for iz in range(ncell):
            filename = root+str(ix)+str(iy)+str(iz)+'.txt'
            if os.path.exists(filename):
                os.remove(filename)

# Create a README file with the header for all the files
filename = outdir+'README.md' ; print filename
with open(filename, 'w') as outf:
    outf.write('# fof_halo_mean_x (Mpc/h), fof_halo_mean_y (Mpc/h), fof_halo_mean_z (Mpc/h), fof_halo_mean_vx (km/s), fof_halo_mean_vy (km/s), fof_halo_mean_vz (km/s), log10(fof_halo_mass /Msun/h), fof_halo_tag \n')
    outf.closed
