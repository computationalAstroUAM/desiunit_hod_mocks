#! /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1 
#PBS -q sciama1.q

code2run=/users/gonzalev/eboss/mock_construction/produce_ascii_files/boxes/fof_boxes.py
#code2run=/users/gonzalev/eboss/mock_construction/produce_ascii_files/boxes/check_fofboxes.py
#code2run=/users/gonzalev/eboss/mock_construction/produce_ascii_files/boxes/mass_bins.py

echo '~~~~~~~~~~~~Run'
echo $code2run         
echo '~~~~~~~~~~~~'

# Run
python $code2run
