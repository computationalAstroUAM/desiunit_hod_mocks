#! /bin/tcsh -f
cd ${PBS_O_WORKDIR}

set snap = ${snap}
set ivol = ${ivol}
set exec = snap2ascii.py

echo '~~~~~~~~~~~~Run'
echo $snap $ivol
echo $exec         
echo '~~~~~~~~~~~~'

# Run
python $exec $snap $ivol
