# The OuterRim simulation  

The OuterRim simulation was run on a cubic volume of side 3000 Mpc$/h$, with dark matter particles with masses of $1.9\cdot 10^9 M_{\odot}/h$. The simulation has 100 snapshots, with 34 between $1<z<3.5$. 

## Friend-of-friends halo information

Between $z=0$ and $z=10$ a friend-of-friend halo catalogue has been contructed assuming $b=0.168$. There are a minimun of 20 particles per halo.

The OuterRim Friend-of-friend halo information has been stored in 110 files ending "#0...#109".

The information in the files "\*fofproperties\*" under each "STEP#" directory is as follow:

```
[data type] Variable name

[i 32] fof_halo_count
[i 64] fof_halo_tag
[f 32] fof_halo_mass
[f 32] fof_halo_center_x
[f 32] fof_halo_center_y
[f 32] fof_halo_center_z
[f 32] fof_halo_mean_x
[f 32] fof_halo_mean_y
[f 32] fof_halo_mean_z
[f 32] fof_halo_mean_vx
[f 32] fof_halo_mean_vy
[f 32] fof_halo_mean_vz
[f 32] fof_halo_vel_disp
```

### Definitions from Katrin Heitmann:

* These files are just FOF halos, no substructure information. While the subfiles do contain subvolumes, they are not contiguous volumes. So in one subfile you can have several cubes of halos that are distributed throughout the full volume.

* The halo count is the number of particles in a halo, in some sense it is redundant with the halo mass (since you can just take that number and multiply it by the particle mass). 

* The halo tag is mostly there to enable identifying halo properties in other files at the same time step. The number is not consistent through the time steps. 

* The halo mass is measured in Msun/h. 

* The fof_halo_center is measured in comoving Mpc/h and it is the potential minimum. 

* The fof_halo_mean is the position of the center of mass. 

* Velocities are comoving peculiar in km/s. 

* The halo velocity dispersion is something I would ignore for now (that was not tested carefully and we never used it). 

## Particle files information

The OuterRim Friend-of-friend halo information has been stored in 110 files ending "#0...#109".

The information in the files "\*particles\*" under each "STEP#" directory is as follow:

```
[data type] Variable name

[f 32] x
[f 32] y
[f 32] z
[f 32] vx
[f 32] vy
[f 32] vz
[i 64] id
[i 64] fof_halo_tag
```

* x,y,z are measured in comoving Mpc/h 
* vx,vy, vz are comoving peculiar in km/s.
* id goes from 0 to the total number of particles.
* fof_halo_tag is the same tag as in the fof files.

# Constructing mock catalogues from teh OuterRim simulation

This repository contains code relevant for the construction of mocks catalogues using the OuterRim simulation.

* **produce_ascii_files**: Contains programs to read the OuterRim simulation files in the genericIO format (both halo and particle information) and produce ASCII files. In order to use this code, the python module generio is required. This can be get following [this directions](http://trac.alcf.anl.gov/projects/genericio).
