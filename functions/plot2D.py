



#---- Modules ----#

import matplotlib.pyplot as plt
import numpy as np



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