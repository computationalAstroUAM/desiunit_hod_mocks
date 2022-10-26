#!/bin/bash

#Lbox = ...                             #If UNIT box: 1000
#Om = ...                               #If UNITSIM cosmology: Om = 0.3089
#z1  = float(sys.argv[3])               #If ELG eBOSS: z = 0.865
#muAcAs = float(sys.argv[4])            #Example: mu10.841_Ac0.0062_As0.01418
#beta = float(sys.argv[5])              #We consider -2 NI  ; -1-0 B ; 0 P ; 0-... NB


PATH1='../DESI_outputs/output_V1/21particles_fsat_beta'
echo "Leer archivos..."
find "$PATH1" -mindepth 0 -type f -name "galaxies_1000Mpc_V1.4more_NFW_*_vfact1.00_beta*_K1.00_vt0pm0_BVG_product_nosubhalos_21particles.dat" $FIRST $SECOND | while read -r file
do
    python3 nersc_tables.py 1000 0.3089 0.865 $FIRST $SECOND
done



