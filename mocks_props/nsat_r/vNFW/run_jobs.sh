#!/bin/bash

Testing=false
sciama=true

logpath=/mnt/lustre/gonzalev/Junk

path2code=/users/gonzalev/eboss/mock_construction/mocks_props/nsat_r/vNFW/

code2run=nsat_r_vNFW.py

K=(1.00 0.40 0.60 0.80)

if $Testing ; then
    K=(0.40)
fi

for ik in ${K[@]} ; do
    name=nsat_K${ik}
    logname=${logpath}/${name}.%j.log

    if ! $sciama ;then
	time python $code2run $ik
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
python $code2run $ik" > $job_file
	
	sbatch $job_file
	rm $job_file
    fi
done

echo 'End of the script'
