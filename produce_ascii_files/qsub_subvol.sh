#! /bin/bash
#PBS -l walltime=00:30:00
#PBS -l nodes=1

code2run=subvol2ascii.py
lcut=200

steps=(203 266 300) #OuterRim: (203 266 300)

logpath=/mnt/lustre/$(whoami)/Junk/$code2run

for step in ${steps[@]} ; do
    logfile=$logpath$step
    rm -f $logfile

    qsub -q sciama1.q -o $logfile -j oe run_subvol.sh -v step=$step,lcut=$lcut,code2run=$code2run

    # Testing
    #qsub -I run_subvol.sh -v step=$step,lcut=$lcut,code2run=$code2run
done

echo 'End of the script'
