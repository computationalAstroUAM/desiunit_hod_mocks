#!/bin/bash

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
