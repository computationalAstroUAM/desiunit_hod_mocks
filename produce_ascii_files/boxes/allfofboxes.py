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

# Minimum values
x0, y0, z0 = [np.zeros(ncell) for i in range(3)]
for i in range(3):
    x0[i] = i*cell
    y0[i] = i*cell
    z0[i] = i*cell

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
    outf.write('# fof_halo_center_x (Mpc/h), fof_halo_center_y (Mpc/h), fof_halo_center_z (Mpc/h), log10(fof_halo_mass /Msun/h), fof_halo_tag \n')
    outf.closed

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
    
    size = len(xc)
    for i,x in enumerate(xc):
        y = yc[i] ; z = zc[i] 
        ix = int(x*float(ncell)/lbox) #; print (ix, x0[ix])
        iy = int(y*float(ncell)/lbox) #; print (iy, y0[iy])
        iz = int(z*float(ncell)/lbox) #; print (iz, z0[iz])
        ioc = ix + ncell*iy +ncell2*iz
        outfile = root+str(ix)+str(iy)+str(iz)+'.txt'

        if (ioc>ncell3 or ioc<0):
            print('ERROR: Cell ',ioc,' of ',ncell3)
            print('ERROR: Cell ',ix,iy,iz)
            print('ERROR: Position ',x,y,z)
            sys.exit()

        if (mass[i]>0.):
            lmass = np.log10(mass[i])
            
            # Shift values so all the boxes start at (0,0,0)
            xs = x - x0[ix]
            ys = y - y0[iy]
            zs = z - z0[iz]

            tofile = ' '.join(map(str, (xs,ys,zs,vx[i],vy[i],vz[i],lmass,tag[i])))

            if os.path.exists(filename):
                with open(outfile, 'a') as outf:
                    outf.write(tofile+' \n') #; print outfile
                    outf.closed
            else:
                with open(outfile, 'w') as outf:
                    outf.write(tofile+' \n') #; print outfile
                    outf.closed

#        # Testing -------------
#        if i>4:
#            break
#    ifile += 1
#    if ifile>4:
#        break
#    #-------------------------

print ('Output in ',root)
