#!/bin/bash

Testing=true
sciama=true

logpath=/mnt/lustre/gonzalev/Junk

path2code=/users/gonzalev/eboss/mock_construction/mocks_props/dprop/

code2run=test.py #dprop.py

typemock=NFW #NFW; part

xx=(0 1 2)
yy=(0 1 2)
zz=(0 1 2)

if $Testing ;then
    #code2run=test.py
    typemock=testNFW
    xx=(0)
    yy=(1)
    zz=(2)
fi

for ix in ${xx[@]} ; do
    for iy in ${yy[@]} ; do
	for iz in ${zz[@]} ; do
	    name=dprop_${ix}${iy}${iz}
	    logname=${logpath}/${name}.%j.log

	    while IFS= read -r mock; do
		if ! $sciama ;then
		    time python $code2run $ix $iy $iz $typemock $mock
		else
		    job_file=${logpath}/${name}.job
		    echo "#!/bin/bash
#
#SBATCH --nodes=1  
#SBATCH --ntasks=1
#SBATCH --time=0-9:00:00
#SBATCH -p sciama4-12.q
#SBATCH --job-name=${name}
#SBATCH -o ${logname}  
#SBATCH -D ${path2code}
#
python $code2run $ix $iy $iz $typemock $mock" > $job_file

		    sbatch $job_file
		    rm $job_file
		    sleep 5s
		fi
	    done < ${typemock}_mocks.txt
	done
    done
done

echo 'End of the script'
