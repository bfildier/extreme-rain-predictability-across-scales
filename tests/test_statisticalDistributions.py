import sys,os
import time

## Own functions
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,os.path.join(os.path.dirname(currentpath),'functions'))

from statisticalDistributions import *

if __name__ == "__main__":

	print("--Test getInvLogRanks")
	n_pts = 1e4
	n_pts_per_bin = 4
	n_bins_per_decade = 10
	start = time.time()
	ranks_invlog = getInvLogRanks(n_pts,n_pts_per_bin,n_bins_per_decade)
	end = time.time()
	print("q's:",ranks_invlog)
	print("time elapsed:", end-start)
	print()

	print("--Test getSample1D")
	x = np.random.random(size=(100,100,100))
	subset_3D = np.random.random(size=(100,100,100))
	subset_3D = np.around(subset_3D)
	subset_3D = subset_3D.astype(bool)
	subset_2D = subset_3D[0,...]
	print('- with 3D slice')
	print("# valid points =",subset_3D.sum())
	sample_1D = getSample1D(x,subset_3D)
	print('sample shape:',sample_1D.shape)
	print('- with 2D slice')
	print("# valid points =",subset_2D.sum()*subset_3D.shape[0])
	sample_1D = getSample1D(x,subset_2D)
	print('sample shape:',sample_1D.shape)
	print()

	print("--Test computePercentilesAndBinsFromRanks")
	x = np.random.random(size=n_pts)
	start = time.time()
	ranks_invlog, centers_invlog, breaks_invlog = computePercentilesAndBinsFromRanks(x,ranks_invlog)
	end = time.time()
	print("x_q's:",centers_invlog)
	print("time elapsed:", end-start)
	print()

	print("--Test computePercentilesAndBinsFromRanks_da")
	x_da = da.from_array(x,chunks = 2)
	start = time.time()
	ranks_invlog, centers_invlog, breaks_invlog = \
		computePercentilesAndBinsFromRanks_da(x,ranks_invlog)
	end = time.time()
	print("x_q's:",centers_invlog)
	print("time elapsed:", end-start)
	print()

	print("--Test defineInvLogPercentiles")
	x = np.random.random(size=n_pts)


	sys.exit(0)
