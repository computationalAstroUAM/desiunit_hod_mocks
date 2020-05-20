#!/bin/bash

sciama=true

logpath=/mnt/lustre/gonzalev/Junk

path2code=/users/gonzalev/eboss/mock_construction/mocks_props/vdis/

code2run=r_vdis.py
#code2run=theta_vdis.py
#code2run=phi_vdis.py

name=vdis
logname=${logpath}/${name}.%j.log

if ! $sciama ;then
    time python $code2run
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
python $code2run" > $job_file
	
    sbatch $job_file
    rm $job_file
fi

echo 'End of the script'
