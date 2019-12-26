#!/bin/bash

logpath=/mnt/lustre/gonzalev/Junk

path2code=/users/gonzalev/eboss/mock_construction/hmf/
code2run=mf.py

steps=(266)
#steps=(300 279 253 241)

for step in ${steps[@]} ; do
    name=hmf_${step}
    logname=${logpath}/${name}.%j.log
    job_file=${logpath}/${name}.job

echo "#!/bin/bash
#
#SBATCH --nodes=1  
#SBATCH --ntasks=1
#SBATCH --time=0-9:00:00
#SBATCH -p himem.q
#SBATCH --job-name=${name}
#SBATCH -o ${logname}  
#SBATCH -D ${path2code}
#
python $code2run $step" > $job_file

    sbatch $job_file
    sleep 5s
    rm $job_file
done

echo 'End of the script'
