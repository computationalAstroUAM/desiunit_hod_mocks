#!/bin/bash

sciama=true

logpath=/mnt/lustre/gonzalev/Junk

path2code=/users/gonzalev/eboss/mock_construction/mocks_props/prob_sat/

code2run=prob_sat.py

name=psat
logname=${logpath}/${name}.%j.log

if ! $sciama ;then
    time python $code2run
else
    job_file=${logpath}/${name}.job
    echo "#!/bin/bash
#
#SBATCH --nodes=1  
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH -p sciama4.q
#SBATCH --job-name=${name}
#SBATCH -o ${logname}  
#SBATCH -D ${path2code}
#
python $code2run" > $job_file
	
    sbatch $job_file
    rm $job_file
fi

echo 'End of the script'
