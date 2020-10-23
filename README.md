# Constructing DESI mock catalogues from the UNIT simulation

This repository contains code relevant for the construction of DESI mocks catalogues using the UNIT simulation.

* **hmf**: Contains code to generate the halo mass function from the OuterRim simulation. [to be adapted to UNIT]

* **biasf**: Code to generate the bias function for the OuterRim simulation. [to be adapted to UNIT]

* **compute_nz**

* **hod_mocks**: Contains programs to find the best HOD values and to populate the simulation accordingly. Code from Martin White can be found here: https://github.com/martinjameswhite/MockingDESI.git

* **mocks_props**: Calculations and plots for exploring the characteristics of different mocks, Nsat vs r, etc.

* **tools**: Contains programs that are useful in all projects, such an example on reading 
in python a GenericIO file.



# DESI specifications
Specifications following the mocks tabulated [here](https://desi.lbl.gov/trac/wiki/Clustering/MockChallenge/post-recon-BAO/stage2).

* UNITSIM-ELG mocks at redshift = 0.9873 (snap 97), densities={25e-4,20e-4,5.5e-4}/(Mpc/h)^-3, bias=1.4