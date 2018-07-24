#! /bin/tcsh -f

set step = (203 266 300)
set stepdir = ('OuterRim_STEP203_z1.433' 'OuterRim_STEP266_z0.865' 'OuterRim_STEP300_z0.656')
set asciidir = '/mnt/lustre/eboss/OuterRim/ascii/'

set nvolmax = 109

@ j = 0
foreach sd ($stepdir)
    @ j = $j + 1
    @ i = 0
    while ($i <= $nvolmax)    
    	set file = ${asciidir}${sd}'/OuterRim_STEP'${step[$j]}'_fofproperties'${i}'.txt'
	if ( !(-f $file) ) then
	    echo 'Not found: '$file
	endif
    	@ i = $i + 1
    end
end

echo 'The end'
