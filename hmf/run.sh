#! /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1 
#PBS -q sciama1.q

root=/users/gonzalev/eboss/mock_construction/hmf/

code2run=${root}mf.py

echo '~~~~~~~~~~~~Run'
echo $code2run         
echo '~~~~~~~~~~~~'

# Run
python $code2run
