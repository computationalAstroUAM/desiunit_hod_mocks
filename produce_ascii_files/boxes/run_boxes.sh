#!/bin/bash

ix=${ix}
iy=${iy}
iz=${iz}
code2run=${code2run}

echo '~~~~~~~~~~~~Run'
echo $code2run   
echo $ix $iy $iz      
echo '~~~~~~~~~~~~'

# Run
python $code2run $ix $iy $iz
