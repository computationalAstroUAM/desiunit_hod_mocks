#! /bin/bash
cd ${PBS_O_WORKDIR}

step=${step}
lcut=${lcut}
code2run=${code2run}

echo '~~~~~~~~~~~~Run'
echo $step $lcut
echo $code2run         
echo '~~~~~~~~~~~~'

# Run
python $code2run $step $lcut
