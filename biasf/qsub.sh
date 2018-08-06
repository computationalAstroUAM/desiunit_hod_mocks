#! /bin/bash
#PBS -l walltime=99:00:00
#PBS nodes=1:ppn=1 
#PBS -q sciama1.q

code2run=/users/gonzalev/eboss/mock_construction/biasf/write_input.py

space=rspace #zspace

xx=(0 1 2)
yy=(0 1 2)
zz=(0 1 2)

for ix in ${xx[@]} ; do
    for iy in ${yy[@]} ; do
	for iz in ${zz[@]} ; do
            qsub run.sh -v space=$space,ix=$ix,iy=$iy,iz=$iz,code2run=$code2run

            ## Testing
            #qsub -I run.sh -v space=$space,ix=$ix,iy=$iy,iz=$iz,code2run=$code2run

	    sleep 5s
	done
    done
done

echo 'End of the script'
