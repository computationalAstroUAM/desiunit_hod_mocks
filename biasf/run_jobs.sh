#!/bin/bash

logpath=/mnt/lustre/gonzalev/Junk

path2code=/users/gonzalev/eboss/mock_construction/biasf/

code2run=write_input.py

space=rspace #zspace

xx=(0) #(0 1 2)
yy=(0) #(0 1 2)
zz=(0) #(0 1 2)

for ix in ${xx[@]} ; do
    for iy in ${yy[@]} ; do
	for iz in ${zz[@]} ; do
	    name=bias_${ix}${iy}${iz}
	    logname=${logpath}/${name}.%j.log

	    job_file=${logpath}/${name}.job
echo "#!/bin/bash
#
#SBATCH --nodes=1  
#SBATCH --ntasks=1
#SBATCH --time=0-9:00:00
#SBATCH -p sciama4.q
#SBATCH --job-name=${name}
#SBATCH -o ${logname}  
#SBATCH -D ${path2code}
#
python $code2run $ix $iy $iz $space" > $job_file

            sbatch $job_file

	    rm $job_file

	    #sleep 5s
	done
    done
done

echo 'End of the script'
