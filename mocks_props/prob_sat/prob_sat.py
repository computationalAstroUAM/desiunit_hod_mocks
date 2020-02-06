import sys,os
import numpy as np
from iotools import check_file

Testing = False

istep = 266

typemock=['NFW','part']

xboxs = ['0','1','2']
yboxs = ['0','1','2']
zboxs = ['0','1','2']

if (Testing):
	typemock=['part']
	xboxs = ['0'] ; yboxs = ['1'] ; zboxs =['2']

# Output path
outpath = '/mnt/lustre/eboss/OuterRim/mocks/nsat_'

# Loop over mocks and boxes
maxn = -999.
for itm, tmock in enumerate(typemock):
	with open(tmock+'_mocks.txt', 'r') as ff:
		mocks = [line.strip() for line in ff]		

	for im, mock in enumerate(mocks):
		print('Mock: {}'.format(mock)) 

		# Loop over boxes
		for xbox in xboxs:
			for ybox in yboxs:
				for zbox in zboxs:
					ibox = xbox+ybox+zbox #; print('Box: {}'.format(ibox))

					# Change the mock names to the box we are working with
					mockfile = mock.replace('mock000','mock'+ibox)
					check_file(mockfile) #; print('Mockfile: {}'.format(mockfile))

					# Output file
					imock = mockfile.split('/')[-1]
					outfile = outpath+tmock+'/'+imock.replace('galaxies','nsat')

					# Read mock catalogue (with or without d*
					# x 0, y 1,z 2(Mpc/h), vx 3, vy 4, vz 5(comoving peculiar in km/s), 
					# lmass 6(np.log10(fof_halo_mass)), cen 7(-1=sat, 0=cen), 
					# dvx 8, dvy 9, dvz 10, dx 11, dy 12, dz 13, tag 14(fof_halo_tag)  
					with open(mockfile, 'rb') as ff: 
						first_line = ff.readline()
						if (len(first_line.strip().split()) == 9):
							itag = 8
						elif (len(first_line.strip().split()) == 15):
							itag = 14
						else:
							print('STOP: unknown file set-up') ; sys.exit()

					ltags = [] ; lstag = []
					with open(mockfile, 'rb') as ff:
						for line in ff:
							thistag = int(line.strip().split()[itag])
							ltags.append(thistag)

							cen = int(line.strip().split()[7])
							if (cen<0):
								lstag.append(thistag)

					# Check that there are more tags than sat. tags
					if (len(ltags)<len(lstag)):
						print('STOP: number tags={}, n sat.={}'.format(
							len(ltags),len(lstag)))
						sys.exit()

					tags = np.sort(np.unique(np.asarray(ltags,dtype=int)))
					stags = np.sort(np.asarray(lstag,dtype=int))
					ltags = [] ; lstag = []
					nsat = np.zeros(shape=(len(tags)))

					iis = 0  ; continue_loop = True
					for ii,tag in enumerate(tags):
						if (not continue_loop): break

						if (stags[iis] > tag): 
							# No satellite: continue loop
							continue						
						elif (stags[iis] == tag):
							# Satellites: loop
							nsat[ii] = nsat[ii] + 1
							iis += 1
							if (iis>len(stags)-1):
								continue_loop = False
								break

							while (stags[iis] == tag):
								nsat[ii] = nsat[ii] + 1
								iis += 1
								if (iis>len(stags)-1):
									continue_loop = False
									break
						else:
							sys.exit('STOP: Not properly coded the finding nsat!')
							

					# Write Output
					tofile = zip(tags,nsat)
					with open(outfile, 'w') as outf:
						np.savetxt(outf,tofile,fmt=('%i %i'))    
						outf.closed
					#print('Output w tag,nsat: {}'.format(outfile))
	print('Outputpath w tag,nsat: {}'.format(outpath+tmock+'/'))
