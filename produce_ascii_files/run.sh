#! /bin/bash
cd ${PBS_O_WORKDIR}

step=${step}
ifile=${ifile}
code2run=${code2run}

echo '~~~~~~~~~~~~Run'
echo $step $ifile 
echo $code2run         
echo '~~~~~~~~~~~~'

# Run
python $code2run $step $ifile