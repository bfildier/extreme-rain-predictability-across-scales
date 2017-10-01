import sys,os
import time

## Own functions
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,os.path.join(os.path.dirname(currentpath),'functions'))

from statisticalDistributions import *

if __name__ == "__main__":

	print("--Test getInvLogRanks")
	n_pts = 1e7
	n_pts_per_bin = 4
	n_bins_per_decade = 10
	start = time.time()
	max_rank, ranks_IL = getInvLogRanks(n_pts,n_pts_per_bin,n_bins_per_decade)
	end = time.time()
	print("q's:",ranks_IL)
	print("q max:",max_rank)
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
	print("* With a numpy array")
	x = np.random.random(size=n_pts)
	print("array is np.ndarray:",isinstance(x,np.ndarray))
	y = x.copy()
	start = time.time()
	centers_invlog, breaks_invlog = computePercentilesAndBinsFromRanks(x,ranks_IL)
	end = time.time()
	print("x_q's:",centers_invlog)
	print("bin edges:",breaks_invlog)
	print("time elapsed:", end-start)
	print()

	print("* With a dask array")
	y_da = da.from_array(y,chunks = n_pts/4)
	print("array is da.core.Array:",isinstance(y_da,da.core.Array))
	start = time.time()
	centers_invlog_da, breaks_invlog_da = \
		computePercentilesAndBinsFromRanks(y_da,ranks_IL)
	end = time.time()
	print("x_q's:",centers_invlog_da)
	print("bin edges:",breaks_invlog_da)
	print("time elapsed:", end-start)
	print()
	print("* Relative errors between these two versions:")
	print(np.divide(np.array(centers_invlog_da)-centers_invlog,centers_invlog))
	print()

	print("--Test definePercentilesOnInvLogQ")
	print("* With a numpy array")
	start = time.time()
	max_rank, ranks_invlog, centers_invlog, breaks_invlog = \
		definePercentilesOnInvLogQ(x)
	end = time.time()
	print("time elapsed:", end-start)
	print("* With a dask array")
	start = time.time()
	max_rank, ranks_invlog, centers_invlog, breaks_invlog = \
		definePercentilesOnInvLogQ(y_da)
	end = time.time()
	print("time elapsed:", end-start)
	print()

	print("--Test defineLogBins")
	sample = np.random.random(1e5)
	print("Corresponding bins edges:",
		defineLogBins(sample)[1])
	print()

	print("--Test defineLinearBins")
	sample = np.random.random(1e5)
	print("Corresponding bins edges:",
		defineLinearBins(sample)[1])
	print()

	print("--Test computePercentileRanksFromBins")
	centers = np.arange(0,1,.1)
	print("centers =",centers)
	sample = np.random.random(1e5)
	print("Corresponding percentile ranks:",
		computePercentileRanksFromBins(sample,centers))
	print()


	sys.exit(0)
