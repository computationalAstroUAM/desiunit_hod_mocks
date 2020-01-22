#!/bin/bash

Testing=true

logpath=/mnt/lustre/gonzalev/Junk

path2code=/users/gonzalev/eboss/mock_construction/mocks_props/

code2run=dprop.py

xx=(0 1 2)
yy=(0 1 2)
zz=(0 1 2)

if $Testing ;then
    xx=(0)
    yy=(1)
    zz=(2)
fi

for ix in ${xx[@]} ; do
    for iy in ${yy[@]} ; do
	for iz in ${zz[@]} ; do
	    name=bias_${ix}${iy}${iz}
	    logname=${logpath}/${name}.%j.log

	    if $Testing ;then
		python $code2run $ix $iy $iz
	    else
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
python $code2run $ix $iy $iz" > $job_file

                sbatch $job_file

		rm $job_file

		sleep 5s
	    fi
	done
    done
done

echo 'End of the script'
