import sys,os
import numpy as np
from iotools import check_file
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec
#from distinct_colours import get_distinct
#import mpl_style
#plt.style.use(mpl_style.style1)

Testing = False

istep = 266
redshift = 0.865

typemock=['NFW','part']
lsty=['-',':']

xboxs = ['0','1','2']
yboxs = ['0','1','2']
zboxs = ['0','1','2']

if (Testing):
	xboxs = ['0'] ; yboxs = ['1'] ; zboxs =['2']

# Bins in distance (Mpc/h)
vmin = -3000. ; vmax = -vmin ; dv = 20.
vbins = np.arange(vmin,vmax,dv)
vhist = vbins +dv*0.5 #; print(len(vhist)) ; sys.exit() 

# Figure 
fig = plt.figure(figsize=(7,4))
xtit = "$v_{r}({\\rm km/s})$"
ytit = "$P_v$"
xmin = -1500. ; xmax = -xmin
ymin = -0.0001 ; ymax = 0.0107
ddy = 0.08*(ymax-ymin) 

ax = fig.add_subplot(111)
ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax)
ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)

## Count mocks to set colour
#colmax = 0 
#for tmock in typemock:
#	count = 0
#	with open(tmock+'_mocks.txt','r') as ff:
#		for line in ff: count +=1
#	if(count>colmax): colmax = count
#if((colmax % 2)!= 0):
#	print('STOP: odd number of catalogues')
#	sys.exit()
#else:
#	colmax = colmax/2 #taken into account mocks with vinfall
#cols = get_distinct(colmax)

# Santi's colours                                                                     
cols = ['red','green','blue','#bcbd22'] 

# Loop over mocks and boxes
minv = 999. ; maxv = -999.
plot_lines = [] ; plot_colors = [] ; inlegs = [] ; alphas = []
for itm, tmock in enumerate(typemock):
	# Initialise the matrix to store the number of sat at different distances
	with open(tmock+'_mocks.txt', 'r') as ff:
		mocks = [line.strip() for line in ff]

	icol = -1
	for im, mock in enumerate(mocks):
		avgv = 0. ; ntot = 0.

		lsty0 = lsty[itm]
		inleg = tmock 

		vfact1 = mock.split('vfact')[1]
		vfact = '$\\alpha_{v}=$'+vfact1.split('_')[0]

		vt1 = mock.split('vt')[1]
		vt = vt1.split('pm')[0]
		if (float(vt)>0.):
			inleg = tmock+', $v^{\\rm infall}$' 
			if (itm == 0):
				lsty0 = '--' 
			else:
				lsty0 = '-.'
		else:
			icol += 1 #Assumes the names go: vt0, vt500

		nsatv = np.zeros(shape=(len(vbins)))
		for xbox in xboxs:
			for ybox in yboxs:
				for zbox in zboxs:
					ibox = xbox+ybox+zbox #; print('Box: {}'.format(ibox))

					# Change the mock names to the box we are working with
					mockfile = mock.replace('mock000','mock'+ibox)
					check_file(mockfile) ; print('Mockfile: {}'.format(mockfile))

					# Read mock catalogue
					# x 0, y 1,z 2(Mpc/h), vx 3, vy 4, vz 5(comoving peculiar in km/s), 
					# lmass 6(np.log10(fof_halo_mass)), cen 7(-1=sat, 0=cen), 
					# dvx 8, dvy 9, dvz 10, dx 11, dy 12, dz 13, tag 14(fof_halo_tag)  
					ldx = [] ; ldy = [] ; ldz = []
					ldvx = [] ; ldvy = [] ; ldvz = []
					with open(mockfile, 'rb') as ff:
						for line in ff:
							cen = int(line.strip().split()[7])
							if (cen<0):
								ldvx.append(float(line.strip().split()[8]))
								ldvy.append(float(line.strip().split()[9]))
								ldvz.append(float(line.strip().split()[10]))
								ldx.append(float(line.strip().split()[11]))
								ldy.append(float(line.strip().split()[12]))
								ldz.append(float(line.strip().split()[13]))
					dvx = np.asarray(ldvx,dtype=float)/(1+redshift)
					dvy = np.asarray(ldvy,dtype=float)/(1+redshift)
					dvz = np.asarray(ldvz,dtype=float)/(1+redshift)
					dx = np.asarray(ldx,dtype=float)
					dy = np.asarray(ldy,dtype=float)
					dz = np.asarray(ldz,dtype=float)
					ldx = [] ; ldy = [] ; ldz = []
					ldvx = [] ; ldvy = [] ; ldvz = []

					# Get vr for satellite galaxies
					vsat = (dx*dvx + dy*dvy + dz*dvz)/np.sqrt(dx*dx + dy*dy + dz*dz)
					if(np.min(vsat)<minv) : minv = np.min(vsat)
					if(np.max(vsat)>maxv) : maxv = np.max(vsat)

					# Average v2
					avgv = avgv + np.sum(vsat*vsat) 
					ntot = ntot + len(vsat)

					# Histogram
					H, bin_edges = np.histogram(vsat, bins=np.append(vbins,vmax))
					nsatv = nsatv + H
					print('    Nmin={:.2f}, Nmax={:.2f}'.format(np.min(nsatv),np.max(nsatv)))

		area = np.sum(nsatv)*dv
		if (itm==0):
			col0=cols[icol]
		else:
			col0=cols[1]

		l1, = ax.plot(vhist,(nsatv/area),label=inleg,
					  linestyle=lsty0,color=col0)

		if (itm==0 and vt=='0'):
			plot_colors.append(l1)
			alphas.append(vfact)
		if (im<2):
			plot_lines.append(l1)
			inlegs.append(inleg)

		# Arrows to show sqrt(<v^2>)
		avgv = np.sqrt(avgv/ntot)
		ax.arrow(avgv,ymax-ddy,0,ddy,color=col0,ls=lsty0,             
                 length_includes_head=True,                                           
                 head_width=ddy*0.3, head_length=ddy*0.3)  

print('    Vmin={:.2f}, Vmax={:.2f} km/s'.format(minv,maxv))

# Legend
legend1 = ax.legend(plot_lines, inlegs, loc=2)
legend1.draw_frame(False)
for h in legend1.legendHandles:
    h.set_color('k')

plt.gca().add_artist(legend1)

leg = ax.legend(plot_colors, alphas, loc=1,
				handlelength=0, handletextpad=0)
leg.draw_frame(False)
for ii,text in enumerate(leg.get_texts()):
    text.set_color(cols[ii])
for item in leg.legendHandles:
    item.set_visible(False)

# Save figure
plotfile = '/mnt/lustre/eboss/OuterRim/mocks/plots/vdis/r_vdis_'+str(istep)+'.png'
fig.savefig(plotfile,dpi=300)
print('Output: ',plotfile)
