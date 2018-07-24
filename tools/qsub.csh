#! /bin/tcsh -f
#PBS -l walltime=24:00:00
#PBS -l nodes=1

set logpath = /mnt/lustre/$user/Junk/ascii

qsub -q sciama1.q -o $logpath -j oe run.csh 
