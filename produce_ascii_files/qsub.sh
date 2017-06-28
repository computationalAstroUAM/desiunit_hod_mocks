#! /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1

steps=(203 266 300) #OuterRim: (203 266 300)

#code2run=step2ascii.py
#nfiles=110 # OuterRim has a total of 110 fof files: 0-109

code2run=particles2ascii.py
nfiles=256 # OuterRim has a total of 256 fof files: 0-255

logpath=/mnt/lustre/$(whoami)/Junk/$code2run

for step in ${steps[@]} ; do
    for ((ifile=0;ifile<$nfiles;ifile++)) ; do
	logfile=$logpath$step.$ifile
	rm -f $logfile

	qsub -q sciama1.q -o $logfile -j oe run.sh -v step=$step,ifile=$ifile,code2run=$code2run

	# Testing
	#qsub -I run.sh -v step=$step,ifile=$ifile,code2run=$code2run
    done
done

echo 'End of the script'
