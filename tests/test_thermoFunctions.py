import sys,os
import numpy as np

## Own functions
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,os.path.join(os.path.dirname(currentpath),'functions'))

from thermoFunctions import *

if __name__ == "__main__":

	print("--Test airDensity:")
	temp = np.linspace(300,200,num=50)
	pres = np.linspace(10000,300,num=50)
	shum = np.linspace(0.04,0.01,num=50)
	print(airDensity(temp,pres,shum))
	print()

	print("--Test saturationVaporPressure")
	print(saturationVaporPressure(temp))
	print()

	print("--Test saturationSpecificHumidity")
	print(saturationSpecificHumidity(temp,pres))

	sys.exit(0)
