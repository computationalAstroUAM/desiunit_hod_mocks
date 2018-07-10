#!/bin/bash

step=266 #OuterRim: (203 266 300)
lcut=1000
code2run=/users/gonzalev/eboss/mock_construction/produce_ascii_files/subvol2ascii.py


echo '~~~~~~~~~~~~Run'
echo $step $lcut
echo $code2run         
echo '~~~~~~~~~~~~'

# Run
python $code2run $step $lcut
