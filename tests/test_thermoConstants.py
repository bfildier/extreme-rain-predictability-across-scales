import sys,os
import numpy as np

## Own functions
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,os.path.join(os.path.dirname(currentpath),'functions'))

from thermoConstants import *

if __name__ == "__main__":

	vars = locals().copy()
	n = max(map(len,vars))
	print("Variables in the local environment:")
	for var in vars.keys():

		print(('{:%ds} :'%n).format(var), vars[var])