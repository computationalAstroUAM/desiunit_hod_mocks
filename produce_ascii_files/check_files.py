import sys, os, glob, getpass
import numpy as np

steps = np.array([203, 266, 300])

check_fof = False
check_particle = True

for istep in steps:
    halodir = '/mnt/lustre/eboss/OuterRim/OuterRim_sim/'
    nroot = halodir+'HaloCatalog/STEP'+str(istep)  

    # Get the conversion between the name of the time step and redshift
    step = np.genfromtxt(halodir+'step_redshift.txt',usecols=0,dtype=int)
    redshift = np.genfromtxt(halodir+'step_redshift.txt',usecols=1) 
    zz = redshift[np.where(step == istep)]
    ascii = halodir+'ascii/OuterRim_STEP'+str(istep)+'_z'+str(zz[0])+'/' 

    # FoF
    files = glob.glob(nroot+'/*'+str(istep)+'*.fofproperties#*')
    nfiles = len(files) # Input
    print nfiles,' fof input files, step=',istep
    nums = [int(i.split('#', 1)[1]) for i in files]
    print '      from',min(nums),' to ',max(nums)
    
    if (check_fof):
        files = glob.glob(ascii+'OuterRim_STEP'+str(istep)+'_fofproperties_*')
        if (nfiles != len(files)):
            onums = np.asarray([int(i.split('#', 1)[1]) for i in files])
            onums = np.sort(onums)
            missing = set(range(nfiles)).difference(onums)
            print 'STOP: FoF missing file(s) ',missing
        else:
            onums = np.asarray([int(os.stat(i).st_size) for i in files])
            ind = np.where(onums < 10)
            if (np.shape(ind)[1]>0):
                idx = np.asarray([int(i) for i in ind[0]])
                empty = [files[i] for i in idx]
                print 'STOP: Empty file(s) ',empty
            else:
                print '      GOOD: ascii halo files produced, e.g. '
                print '     ',files[0]

                # Further test
                further= True
                if (further):
                    for infile in files:
                        with open(infile) as ff:
                            fof = []
                            for line in ff: 
                                #count, tag, mass, xc, yc, zc, xm, ym, zm, vx, vy, vz
                                vz = np.append(fof, line.strip().split()[11])
                                if (len(vz)>100.): break
                        print 'OK',infile

    # Particles
    files = glob.glob(nroot+'/*'+str(istep)+'*.haloparticles#*')
    nfiles = len(files) # Input
    print len(files),' particles input files, step=',istep
    nums = [int(i.split('#', 1)[1]) for i in files]
    print '      from',min(nums),' to ',max(nums)

    if (check_particle):
        files = glob.glob(ascii+'OuterRim_STEP'+str(istep)+'_particles_*')
        if (nfiles != len(files)):
            onums = np.asarray([int(i.split('#', 1)[1]) for i in files])
            onums = np.sort(onums)
            missing = set(range(nfiles)).difference(onums)
            print 'STOP: FoF missing file(s) ',missing
        else:
            onums = np.asarray([int(os.stat(i).st_size) for i in files])
            ind = np.where(onums < 10)
            if (np.shape(ind)[1]>0):
                idx = np.asarray([int(i) for i in ind[0]])
                empty = [files[i] for i in idx]
                print 'STOP: Empty file(s) ',empty
            else:
                print '      GOOD: ascii particle files produced, e.g. '
                print '     ',files[0]

                # Further test
                further= True
                if (further):
                    for infile in files:
                        with open(infile) as ff:
                            fof = []
                            for line in ff: #x,y,z,vx,vy,vz,ids,fof
                                fof = np.append(fof, line.strip().split()[7])
                                if (len(fof)>100.): break
                        print 'OK',infile
