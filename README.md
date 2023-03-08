# Constructing eBOSS and DESI mock catalogues from the UNIT simulation

The code is based on [Avila+2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.499.5486A/abstract). Please cite this paper if you use the code.

This repository contains code relevant for the construction of eBOSS and DESI mock catalogues using the UNIT simulation.

* **hmf**: Contains code to generate the halo mass function from the OuterRim simulation. [to be adapted to UNIT]

* **biasf**: Code to generate the bias function for the OuterRim simulation. [to be adapted to UNIT]

* **fit_mocks**: Contains code for obtaining HOD parameters As, Ac, mu, fsat

* **hod_mocks**: Contains programs to find the best HOD values and to populate the simulation accordingly. 

* **mocks_props**: Calculations and plots for exploring the characteristics of different mocks, Nsat vs r, etc.

* **nersc_tables.py**: to generate tables for NERSC 

* **DESI_Paramters**: contains a file indicating the values of several HOD parameters used for the fitting of DESI ELG fuji z = 0.8-1.1 data


# DESI specifications
Specifications following the mocks tabulated [here](https://desi.lbl.gov/trac/wiki/Clustering/MockChallenge/post-recon-BAO/stage2).
* UNITSIM-ELG mocks at redshift = 0.9873 (snap 97), densities={25e-4,20e-4,5.5e-4}/(Mpc/h)^-3, bias=1.4

#eBOSS specifications
* UNITSIM-ELG mocks at redhsift = 0.8594 (snap 100), densities = (2.1e-4, 6·2.1e-4),bias = 1.37

