"""Module plotMaps

"""


#---- Modules ----#
import warnings
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm
import numpy as np
import numpy.ma as ma
from math import log10


#---- Functions ----#

def plotMapTropics(lon2D,lat2D,v,v_min=None,v_max=None,mode=None,above_zero=False,\
                   cmap=plt.cm.Greens,title=None,maskcontinents=False,plotfunction='contourf'):
    
    """Assumes that the data is 2D and ranges 30S-30N."""

    v2plot = v
    if mode == "log" or above_zero:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            v2plot = ma.masked_where(v <= 1e-14,v)
    
    if v_min is None:
        v_min = np.nanmin(v2plot)
    if v_max is None:
        v_max = np.nanmax(v2plot)

    nlev = 11
    dl = nlev // 5

    fig = plt.figure(figsize=(12,2))
    
    map = Basemap(projection='cyl',lat_0=0,lon_0=180,llcrnrlon=0,llcrnrlat=-30,\
        urcrnrlon=357,urcrnrlat=30)
    if mode == "log":
        levs = np.logspace(log10(v_min),log10(v_max),nlev)
        norm = LogNorm()
    else:
        levs = np.linspace(v_min,v_max,nlev)
        norm = None
    if plotfunction == 'contourf':
        map.contourf(lon2D,lat2D,v2plot,levels=levs,cmap=cmap,norm = norm)    
        #         map.contourf(lon2D,lat2D,freqMap,cmap=plt.cm.Greens)
    elif plotfunction == 'pcolormesh':
        map.pcolormesh(lon2D,lat2D,v2plot,cmap=cmap,norm = norm,
            vmin=v_min,vmax=v_max)
    map.drawparallels(range(-90, 100, 30),labels=[1,0,0,1])
    map.drawmeridians(range(0,400,90),labels=[1,0,0,1])
    map.drawcoastlines()
    if maskcontinents:
        map.fillcontinents(color='k')

    plt.title(title)
    plt.colorbar(ticks=levs[::dl],pad=0.02,fraction=0.085)




