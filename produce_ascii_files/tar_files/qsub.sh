#! /bin/bash
#PBS -l walltime=48:00:00
#PBS -l nodes=1

logpath=/mnt/lustre/$(whoami)/Junk/

indir='/mnt/lustre/eboss/OuterRim/OuterRim_sim/ascii/'
outdir='/mnt/astro3/sftp/ebossguest/downloads/ascii/'

#reds=('OuterRim_STEP203_z1.433' 'OuterRim_STEP266_z0.865' 'OuterRim_STEP300_z0.656')
reds=('OuterRim_STEP300_z0.656')
#whats=('fofproperties' 'particles')
whats=('particles')

for red in ${reds[@]} ; do
    for what2tar in ${whats[@]} ; do
	logfile=$logpath${red}_tar_$what2tar

	qsub -q sciama1.q -o $logpath -j oe run.sh -v indir=$indir,outdir=$outdir,red=$red,what2tar=$what2tar

	# Testing
	#qsub -I run.sh -v indir=$indir,outdir=$outdir,red=$red,what2tar=$what2tar 
    done
done
