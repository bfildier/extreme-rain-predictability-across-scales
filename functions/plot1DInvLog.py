"""Module plot1DInvLog

Functions to plot curves on inverse-logarithmic x-axis (extremes).
"""


#---- Modules ----#

import matplotlib.pyplot as plt
import numpy as np
from math import log10,ceil


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

def subplotRanksILog(ax,ranks,y,col=None,ltype=None,alpha=None,transformX=False):
    
    ax.set_xscale('log')

    # define x-axis
    x = np.flipud(1./(1-ranks/100.))
    # plot
    if isinstance(y,list):
        for i in range(len(y)):
            lt = [ltype[i] if ltype is not None else '-'][0]
            a = [alpha[i] if alpha is not None else 1][0]
            c = [col[i] if col is not None else 1][0]
            ax.plot(x,y[i],c=c,alpha=a,linestyle=lt)
    else:
        ax.plot(x,y,c=col,alpha=alpha,linestyle=ltype)

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
    # plot
    x0 = x[iQ_lims[0]]
    x1 = x[iQ_lims[1]]
    ax.axvspan(x0,x1,color = '0.75',alpha=alpha)
    
    # transform x-axis
    if transformX:
        transformXaxisIL(ax,x)

def addZeroLine(ax,x):

    ax_line = ax.twinx()
    subplotRanksILog(ax_line,x,
                     np.zeros(x.size),
                     col='gray',ltype=':',transformX=False)
    ax_line.yaxis.set_ticks_position('none')
    ax_line.yaxis.set_ticklabels('')
    ax_line.set_ylim(ax.get_xlim())

