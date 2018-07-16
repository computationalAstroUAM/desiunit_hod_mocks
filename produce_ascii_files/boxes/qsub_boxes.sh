#! /bin/bash
#PBS -l walltime=72:00:00
#PBS nodes=1:ppn=1 
#PBS -q sciama1.q

code2run=/users/gonzalev/eboss/mock_construction/produce_ascii_files/boxes/fof_box.py

xx=(0) #(0 1 2)
yy=(1) #(0 1 2)
zz=(2) #(0 1 2)

for ix in ${xx[@]} ; do
    for iy in ${yy[@]} ; do
	for iz in ${zz[@]} ; do
            #qsub -q sciama1.q -l walltime=74:00:00 -o $logfile -e $logfile.e run_subvol.sh -v ix=$ix,iy=$iy,iz=$iz,code2run=$code2run

            # Testing
            qsub -I run_boxes.sh -v ix=$ix,iy=$iy,iz=$iz,code2run=$code2run

	    sleep 5s #1h
	done
    done
done

echo 'End of the script'
