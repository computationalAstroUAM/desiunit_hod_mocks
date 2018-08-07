import glob
import os, time

runfiles = glob.glob('/mnt/lustre/gonzalev/Junk/xi_mass_bin/*.sh')

for ff in runfiles:
	os.system('qsub '+ff)
	time.sleep(5)
