#! /bin/bash
#PBS -l walltime=24:00:00
#PBS nodes=1:ppn=1 
#PBS -q sciama1.q

code2run=subvol2ascii.py
lcut=1000

steps=(266) #OuterRim: (203 266 300)

logpath=/mnt/lustre/$(whoami)/Junk/$code2run

for step in ${steps[@]} ; do
    logfile=$logpath$step
    rm -f $logfile

    #qsub -q sciama1.q -l walltime=24:00:00 -o $logfile -e $logfile.e run_subvol.sh -v step=$step,lcut=$lcut,code2run=$code2run

    qsub run_subvol.sh -v step=$step,lcut=$lcut,code2run=$code2run

    # Testing
    #qsub -I run_subvol.sh -v step=$step,lcut=$lcut,code2run=$code2run
done

echo 'End of the script'
