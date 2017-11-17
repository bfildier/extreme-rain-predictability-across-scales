"""Module plot1DInvLog

Functions to plot curves on inverse-logarithmic x-axis (extremes).
"""


#---- Modules ----#

import matplotlib.pyplot as plt
import numpy as np
from math import log10,ceil
from matplotlib.patches import Polygon


#---- Functions ----#

def transformXaxisIL(ax,x):
    
    # reverse x-axis
    ax.invert_xaxis()
    # rename ticks
    labels = [item.get_text() for item in ax.get_xticklabels()]
    n = ceil(log10(x.max()))
    N = len(labels)
    for i in range(1,N):
        labels[-n+i-4] = str(100*(1-10**(-n+i-1)))
        if -n+i-1 == 0:
            break
    ax.set_xticklabels(labels)

def subplotRanksILog(ax,ranks,y,col=None,ltype=None,linewidth=None,alpha=None,transformX=False):
    
    ax.set_xscale('log')

    # define x-axis
    x = np.flipud(1./(1-ranks/100.))
    # plot
    if isinstance(y,list):
        for i in range(len(y)):
            lt = ltype[i] if ltype is not None else '-'
            a = alpha[i] if alpha is not None else 1
            c = col[i] if col is not None else 1
            lw = linewidth[i] if linewidth is not None else 1.5
            ax.plot(x,y[i],c=c,alpha=a,linestyle=lt,linewidth=lw)
    else:
        ax.plot(x,y,c=col,alpha=alpha,linestyle=ltype,linewidth=linewidth)

    # transform x-axis
    if transformX:
        transformXaxisIL(ax,x)
    
def subplotYShadingRanksILog(ax,ranks,y_BCs,col,alpha=0.2,transformX=False):
    
    ax.set_xscale('log')
    
    # define x-axis
    x = np.flipud(1./(1-ranks/100.))
    # plot
    y1 = y_BCs[0]
    y2 = y_BCs[1]
    ax.fill_between(x, y1, y2, where=y2 >= y1, facecolor=col,alpha=alpha, interpolate=True)
    
    # transform x-axis
    if transformX:
        transformXaxisIL(ax,x)

def subplotXShadingRanksILog(ax,ranks,iQ_lims,alpha=0.2,col='0.75',transformX=False):

    ax.set_xscale('log')
    
    # define x-axis
    x = np.flipud(1./(1-ranks/100.))
    if iQ_lims[0] >= x.size:
        return
    x0 = x[iQ_lims[0]]
    if iQ_lims[1] >= x.size:
        x1 = x[-1]
    else:
        x1 = x[iQ_lims[1]]
    # plot
    ax.axvspan(x0,x1,color = '0.75',alpha=alpha)
    
    # transform x-axis
    if transformX:
        transformXaxisIL(ax,x)

def addXHatchRanksILog(ax,ranks,iQ_lims,transformX=False):

    ax.set_xscale('log')

    x = np.flipud(1./(1-ranks/100.))
    ax.add_patch(Polygon([[x[iQ_lims[0]], ax.get_ylim()[0]],\
                          [x[iQ_lims[1]], ax.get_ylim()[0]],\
                          [x[iQ_lims[1]], ax.get_ylim()[1]],\
                          [x[iQ_lims[0]], ax.get_ylim()[1]]],\
                          closed=True, fill=False, hatch='//',linewidth=0,
                          color = 'gray'))
    # transform x-axis
    if transformX:
        transformXaxisIL(ax,x)

def addZeroLine(ax,x):

    ax_line = ax.twinx()
    subplotRanksILog(ax_line,x,
                     np.zeros(x.size),
                     col='gray',ltype='-',linewidth=0.8,transformX=False)
    ax_line.yaxis.set_ticks_position('none')
    ax_line.yaxis.set_ticklabels('')
    ax_line.set_ylim(ax.get_ylim())

