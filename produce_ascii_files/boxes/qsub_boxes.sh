#! /bin/bash
#PBS -l walltime=99:00:00
#PBS nodes=1:ppn=1 
#PBS -q sciama1.q

code2run=/users/gonzalev/eboss/mock_construction/produce_ascii_files/boxes/fof_box.py

xx=(0) #0 1 2)
yy=(2) #0 1 2)
zz=(2) #0 1 2)

for ix in ${xx[@]} ; do
    for iy in ${yy[@]} ; do
	for iz in ${zz[@]} ; do
            qsub run_boxes.sh -v ix=$ix,iy=$iy,iz=$iz,code2run=$code2run

            # Testing
            #qsub -I run_boxes.sh -v ix=$ix,iy=$iy,iz=$iz,code2run=$code2run

	    sleep 90s
	done
    done
done

echo 'End of the script'
