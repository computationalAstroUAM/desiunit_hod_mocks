#!/bin/bash

#Lbox = ...                             #If UNIT box: 1000
#Om = ...                               #If UNITSIM cosmology: Om = 0.3089
#z1  = float(sys.argv[3])               #If ELG eBOSS: z = 0.865
#muAcAs = float(sys.argv[4])            #Example: mu10.841_Ac0.0062_As0.01418
#beta = float(sys.argv[5])              #We consider -2 NI  ; -1-0 B ; 0 P ; 0-... NB




python3 nersc_tables.py 1000 0.3089 0.865 mu11.114_Ac0.0109_As0.02871 -0.2

python3 nersc_tables.py 1000 0.3089 0.865 mu10.841_Ac0.0062_As0.01418 -0.2
python3 nersc_tables.py 1000 0.3089 0.865 mu10.841_Ac0.0062_As0.01418 -0.16
python3 nersc_tables.py 1000 0.3089 0.865 mu10.841_Ac0.0062_As0.01418 -0.12
python3 nersc_tables.py 1000 0.3089 0.865 mu10.841_Ac0.0062_As0.01418 -0.08
python3 nersc_tables.py 1000 0.3089 0.865 mu10.841_Ac0.0062_As0.01418 -0.04
python3 nersc_tables.py 1000 0.3089 0.865 mu10.841_Ac0.0062_As0.01418 0
#python3 nersc_tables.py 1000 0.3089 0.865 mu10.841_Ac0.0062_As0.01418 0.04
#python3 nersc_tables.py 1000 0.3089 0.865 mu10.841_Ac0.0062_As0.01418 0.08
#python3 nersc_tables.py 1000 0.3089 0.865 mu10.841_Ac0.0062_As0.01418 0.12
#python3 nersc_tables.py 1000 0.3089 0.865 mu10.841_Ac0.0062_As0.01418 0.16
#python3 nersc_tables.py 1000 0.3089 0.865 mu10.841_Ac0.0062_As0.01418 0.2




