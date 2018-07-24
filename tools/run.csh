#! /bin/tcsh -f
cd ${PBS_O_WORKDIR}

set exec = example_2ascii.py

# Run
python $exec 


