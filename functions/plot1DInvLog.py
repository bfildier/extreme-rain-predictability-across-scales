"""Module plot1DInvLog

Functions to plot curves on inverse-logarithmic x-axis (extremes).
"""


#---- Modules ----#

import matplotlib.pyplot as plt
import numpy as np
from math import log10,ceil


#---- Functions ----#

def transformXaxisIL(ax,x):
    
    ax.invert_xaxis() # reverse x-axis
    labels = [item.get_text() for item in ax.get_xticklabels()]
    n = ceil(log10(x.max()))
    N = len(labels)
    for i in range(1,N):
    #     labels[i] = str(100*(1-10**(-n+i-1)))
        labels[-n+i-4] = str(100*(1-10**(-n+i-1)))
        if -n+i-1 == 0:
            break
    ax.set_xticklabels(labels)

def subplotLogRanksILog(ax,Q_IL,y,col,transformX=False):
    
    ax.set_xscale('log')
    ax.set_yscale('log')

    # define x-axis
    x = np.flipud(1./(1-Q_IL[1:-1]/100.))
    # plot
    if isinstance(y,list):
        for i in range(len(y)):
            ax.plot(x,y[i],c=col[i])
    else:
        ax.plot(x,y)

    # transform x-axis
    if transformX:
        transformXaxisIL(ax,x)
    
def subplotYShadingLogRanksILog(ax,Q_IL,y_BCs,col,alpha=0.2,transformX=False):
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # define x-axis
    x = np.flipud(1./(1-Q_IL[1:-1]/100.))
    # plot
    y1 = y_BCs[0]
    y2 = y_BCs[1]
    ax.fill_between(x, y1, y2, where=y2 >= y1, facecolor=col,alpha=alpha, interpolate=True)
    
    # transform x-axis
    if transformX:
        transformXaxisIL(ax,x)

def subplotXShadingLogRanksILog(ax,Q_IL,i_Q_lims,alpha=0.2,transformX=False):
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # define x-axis
    x = np.flipud(1./(1-Q_IL[1:-1]/100.))
    # plot
    x0 = x[i_Q_lims[0]]
    x1 = x[i_Q_lims[1]]
    ax.axvspan(x0,x1,color = '0.75',alpha=alpha)
    
    # transform x-axis
    if transformX:
        transformXaxisIL(ax,x)