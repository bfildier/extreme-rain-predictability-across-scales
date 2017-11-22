



#---- Modules ----#

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patchess



#---- Functions ----#

def subplotMultiscaleVar(ax,var,time_strides,cmap='viridis',vmin=None,vmax=None):

	if vmin is None:
		vmin = np.nanmin(var)
	if vmax is None:
		vmax = np.nanmax(var)

	plt.imshow(var,cmap=cmap,vmin=vmin,vmax=vmax)
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

	plt.colorbar()

## Hatch time_strides/resolutions where the reference value (e.g. sample
## size) is below some threshold
def addHatchBelowThreshold(ax,var_ref,threshold):

    """Assumes that the plot has been made with imshow,
    so that we can refer to grid cells from their indices."""

    N_res, N_times = samplesize.shape
    for i,j in np.ndindex((N_times,N_res)):
        if var_ref[j,i] <= threshold:
            ax.add_patch(patches.Rectangle((i-0.5, j-0.5), 1, 1, 
               fill=True, snap=False, linewidth=0,color='white'))
            ax.add_patch(patches.Rectangle((i-0.5, j-0.5), 1, 1, 
               hatch='//', fill=False, snap=False, linewidth=0.1,color='gray'))