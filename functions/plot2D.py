



#---- Modules ----#

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from matplotlib.colors import LogNorm
import warnings

#---- Own modules ----#
from plot1DInvLog import *
from statisticalDistributions import *


#---- Functions ----#

def subplotMultiscaleVar(ax,var,time_strides,cmap='viridis',vmin=None,vmax=None):

	if vmin is None:
		vmin = np.nanmin(var)
	if vmax is None:
		vmax = np.nanmax(var)

	im = plt.imshow(var,cmap=cmap,vmin=vmin,vmax=vmax)
	ax.set_xlabel('Time scale')
	ax.set_ylabel('Resolution')
	ax.invert_yaxis()

	# t_hours = list(map(inHours,time_strides))
	# x_km = list(map(coarseningFactor,resolutions))
	# time_strides_2D, resolutions_2D = np.meshgrid(t_hours,x_km)
	# plt.pcolor(time_strides_2D,resolutions_2D,norm_var_pr_sc_prQ_Q25_Q35)
	# ## Adjust x labels
	# xticklabels = ['1 hr','6 hr','1 day','8 days']
	# xticks = list(map(inHours,xticklabels))
	# # xticks = [exp(xi) for xi in np.convolve(np.log(xticks),(0.5,0.5),mode='valid')]
	# ax.set_xticks(xticks)
	# ax.set_xticklabels(xticklabels)
	# ## Adjust y labels
	# yticklabels = ["%d km"%(240*i) for i in range(1,len(ax.get_yticks()))]
	# ax.set_yticks(np.convolve(ax.get_yticks(),(0.5,0.5),mode='valid'))
	# ax.set_yticklabels(labels=yticklabels)

	# x labels
	xticklabels = time_strides
	xticks = list(range(len(xticklabels)))
	ax.set_xticks(xticks)
	ax.set_xticklabels(xticklabels)
	# y labels
	yticklabels = ["%d km"%(240*i) for i in range(len(ax.get_yticks()))]
	ax.set_yticklabels(labels=yticklabels)

	# Colorbar
	cb = plt.colorbar(im)

	return cb

## Hatch time_strides/resolutions where the reference value (e.g. sample
## size) is below some threshold
def addHatchBelowThreshold(ax,var_ref,threshold):

    """Assumes that the plot has been made with imshow,
    so that we can refer to grid cells from their indices."""

    N_res, N_times = var_ref.shape
    for i,j in np.ndindex((N_times,N_res)):
        if var_ref[j,i] <= threshold:
            ax.add_patch(patches.Rectangle((i-0.5, j-0.5), 1, 1, 
               fill=True, snap=False, linewidth=0,color='white'))
            ax.add_patch(patches.Rectangle((i-0.5, j-0.5), 1, 1, 
               hatch='//', fill=False, snap=False, linewidth=0.1,color='gray'))



def subplot2DRanksILog(ax,ranksX,ranksY,Z,cmap=plt.cm.RdBu_r,alpha=None,
	transformX=False,transformY=False,range_type='sym_to_one',vmin=None,vmax=None,
	Z_mode='linear',ranksZ=None):

	"""ranksZ only used if Z_mode is invlog."""

	warnings.filterwarnings("ignore", category=RuntimeWarning)

	ax.set_xscale('log')
	ax.set_yscale('log')

	# define axes
	x = np.flipud(1./(1-ranksX/100.))
	y = np.flipud(1./(1-ranksY/100.))
	X,Y = np.meshgrid(x,y)

	def getRange(Z,vmin,vmax):
		# compute display range
		if vmin is None or vmax is None:
			if np.all(np.isnan(Z)):
				return 0.1,10
			if range_type == 'sym_to_one':
				expmax = int(np.log10(np.nanmax(Z)))
				vmax = 10**expmax
				vmin = 1/vmax
			elif range_type == 'full':
				vmin = np.nanmin(Z)
				vmax = np.nanmax(Z)
			elif range_type == 'full_positive':
				vmin = np.nanmin(Z[Z>0])
				vmax = np.nanmax(Z)
		return vmin,vmax

	def plot(X,Y,Z,cmap,alpha,vmin,vmax,Z_mode,ranksZ):
		# mask nan values
		m = np.ma.masked_where(np.isnan(Z),Z)
		if Z_mode == 'linear':
			im = ax.pcolor(X,Y,m, cmap=cmap,alpha=alpha,vmin=vmin,vmax=vmax)
		elif Z_mode == 'log':
			im = ax.pcolor(X,Y,m, cmap=cmap,alpha=alpha,
				norm=LogNorm(vmin=vmin, vmax=vmax))
		elif Z_mode == 'invlog':
			# transformed values to plot
			Z_new = np.empty(Z.shape)
			Z_new[:] = np.nan
			for ind in np.ndindex(Z.shape):
				if not np.isnan(Z[ind]):
					i_rank = indexOfRank(Z[ind],ranksZ)
					# print(Z[ind],i_rank,ranksZ[-1-i_rank])
					Z_new[ind] = ranksZ[-1-i_rank]
			m = np.ma.masked_where(np.isnan(Z_new),Z_new)
			# find range of transformed values
			vmin = 1/(1-np.nanmax(ranksZ)/100)
			vmax = 1
			# revert cmap
			cmap_r = '%s_r'%cmap if not cmap.endswith('_r') else cmap[:-2]
			im = ax.pcolor(X,Y,1/(1-m/100), cmap=cmap_r,alpha=alpha,
				norm=LogNorm(vmin=vmin,vmax=vmax))
		return im

	# plot
	if isinstance(Z,list):
		for i in range(len(y)):
			a = alpha[i] if alpha is not None else 1
			cm = cmap[i] if isinstance(cmap,list) else cmap
			vmin,vmax = getRange(Z[i],vmin,vmax)
			im = plot(X,Y,Z,cmap,alpha,vmin,vmax,Z_mode,ranksZ)
	else:
		# warnings.filterwarnings("ignore", category=RuntimeWarning)
		vmin,vmax = getRange(Z,vmin,vmax)
		im = plot(X,Y,Z,cmap,alpha,vmin,vmax,Z_mode,ranksZ)

	# colorbar
	if range_type == 'sym_to_one':
		cb = plt.colorbar(im,ticks=[vmin,1,vmax])
		cb.ax.set_yticklabels([str(vmin), '1', str(vmax)])
	else:
		if Z_mode == 'invlog':
			cb = plt.colorbar(im)
			tklb = cb.ax.get_yticklabels()
			imax = ranksZ.size // 10
			newtklb = []
			for tk in tklb:
				if tk.get_text() != '':
					i = int(tk.get_text().split('^{')[-1].split('}}')[0])
					newtklb.append(str(100*(1-10**(i-imax))))
				else:
					newtklb.append('')
			cb.ax.set_yticklabels(newtklb)
			cb.ax.invert_yaxis()
		else:
			cb = plt.colorbar(im)


	# transform x-axis
	if transformX:
		transformXaxisIL(ax,x,offset=1)
	if transformY:
		transformYaxisIL(ax,y,offset=1)

	return cb

