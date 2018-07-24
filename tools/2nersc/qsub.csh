#! /bin/tcsh -f
#PBS -l walltime=24:00:00
#PBS -l nodes=1

set logpath = /mnt/lustre/$user/Junk/ascii

set snaps =  (203 266 300) #OuterRim: (203 266 300)

set nvolmax = 20 # OuterRim: 109
set vols = ()
@ i = 0 #0
while ($i <= $nvolmax)
    set vols = ($vols $i)
    @ i = $i + 1
end

foreach snap ($snaps)
    foreach ivol ($vols)
	set logfile = $logpath$snap.$ivol
	rm -f $logfile

	qsub -q sciama1.q -o $logfile -j oe run.csh -v snap=$snap,ivol=$ivol
    end
end

