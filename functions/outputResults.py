



#---- Modules ----#

import sys,os
import numpy as np
import dask.array as da
import datetime as dt
import pandas as pd

#---- Own functions ----#
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,os.path.join(os.path.dirname(currentpath),'functions'))
from statisticalDistributions import *

#---- Functions ----#

# Display seconds in the format "%H:%M:%S"
def time2str(seconds):
    return (dt.datetime(year=2000,month=1,day=1)+
            dt.timedelta(seconds=(seconds))).strftime("%H:%M:%S")

## Save array size and time for execution of script in dataframe
def saveTiming(fulltimingfile,column_label,arraysize,time_elapsed_s,reset_value=False):
    
    time_elapsed = str(dt.timedelta(seconds=time_elapsed_s))
    if os.path.isfile(fulltimingfile):
        timing_df = pd.read_csv(fulltimingfile)
        timing_df.set_index('index',inplace=True) 
        if reset_value or column_label not in timing_df.keys():
            timing_df[column_label] = [1,arraysize,time_elapsed]
        else:
            # Get previous values
            ntimes = int(timing_df[column_label]['ntimes'])
            oldarraysize = int(timing_df[column_label]['arraysize'])
            runtime_str = timing_df[column_label]['runtime']
            # print warning if sizes are different
            if arraysize != oldarraysize:
                print("WARNING: replacing old array size %d with new size %d in %s."%\
                      (oldarraysize,arraysize,fulltimingfile))
            # Update runtime
            runtime_s = (dt.datetime.strptime(runtime_str.split('.')[0],"%H:%M:%S")-\
                         dt.datetime(1900,1,1)).total_seconds()+\
                        float(".%s"%runtime_str.split('.')[1])
            newruntime_s = (ntimes*runtime_s + (t1-t0))/(ntimes+1)
            newruntime = str(dt.timedelta(seconds=newruntime_s))
            print(runtime_str,newruntime)
            # save
            timing_df[column_label] = [ntimes+1,arraysize,newruntime]
    else:
        timing_df = pd.DataFrame({'index':['ntimes','arraysize','runtime'],column_label:[1,arraysize,time_elapsed]})
        timing_df.set_index('index',inplace=True)
    timing_df.to_csv(fulltimingfile)
    return timing_df

## Load all results for later plotting

# Extract varid in an array for all time_strides and resolutions
def getTXVarFromResults(varid,results,time_strides,resolutions,iQ_slice,avg_mode='mean'):
    
    """'avg_mode' can be mean or sum"""  

    nres = len(resolutions)
    nts = len(time_strides)
    varidTX = np.empty((nres,nts))
    varidTX[:] = np.nan

    for its,ires in np.ndindex(nts,nres):

        time_stride = time_strides[its]
        resolution = resolutions[ires]
        
        try:
            if avg_mode == 'mean':
                varidTX[ires,its] = np.nanmean(results[time_stride][resolution][varid][iQ_slice])
            elif avg_mode == 'sum':
                varidTX[ires,its] = np.nansum(results[time_stride][resolution][varid][iQ_slice])
        except KeyError:
            continue
        
    return varidTX