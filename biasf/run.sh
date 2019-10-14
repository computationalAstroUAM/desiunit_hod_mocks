#!/bin/bash

module purge                                                               
module load system intel_comp/2019.2 cpython/3.7.1

space=${space}
ix=${ix}
iy=${iy}
iz=${iz}
code2run=${code2run}

echo '~~~~~~~~~~~~Run'
echo $code2run   
echo $ix $iy $iz $space      
echo '~~~~~~~~~~~~'

# Run
python $code2run $ix $iy $iz $space
