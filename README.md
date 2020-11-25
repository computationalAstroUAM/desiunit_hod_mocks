# Constructing DESI mock catalogues from the UNIT simulation

This repository contains code relevant for the construction of DESI mocks catalogues using the UNIT simulation.

* **hmf**: Contains code to generate the halo mass function from the OuterRim simulation. [to be adapted to UNIT]

* **biasf**: Code to generate the bias function for the OuterRim simulation. [to be adapted to UNIT]

* **fit_mocks**: Contains code for obtaining HOD parameters As, Ac, mu, fsat

* **hod_mocks**: Contains programs to find the best HOD values and to populate the simulation accordingly. 

* **mocks_props**: Calculations and plots for exploring the characteristics of different mocks, Nsat vs r, etc.

* **nersc_tables.py**: to generate tables for NERSC 

* **pruebas.py**: any probles needed to test python programs



# DESI specifications
Specifications following the mocks tabulated [here](https://desi.lbl.gov/trac/wiki/Clustering/MockChallenge/post-recon-BAO/stage2).

* UNITSIM-ELG mocks at redshift = 0.9873 (snap 97), densities={25e-4,20e-4,5.5e-4}/(Mpc/h)^-3, bias=1.4
