"""Module defineDaskParameters

Define parameters related to the use of dask arrays that are common across modules
"""

## Default compute action for dask arrays
da_compute_default = False
## Size of dask array chunks
chunks = 100000