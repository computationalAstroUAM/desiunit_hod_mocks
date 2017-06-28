# Producing ASCII files

The OuterRim simulation is stored in the binary format genericIO. The programs in this directory read in the OuterRim original files and write them as ASCII files.

* **check_files.py**: Checks how many input fof and particle files there are for a given snapshot. Checks that the same number of ascii files exist.

* **snap2ascii.py**: Deals with the halo catalogue. It can be run using **qsub.sh** and **run.sh**.

* **particles2ascii.py**: Deals with the particle information. It can be run using **qsub.sh** and **run.sh**

* **subvol2ascii.py**: Produces a subvolume of the OuterRim simulation for testing purposes. It can be run using **qsub_subvol.sh** and **run_subvol.sh**