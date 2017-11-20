#!/bin/bash

inputipynb=$1
outputname=$2

# Load python environment
#source activate py36
module load python/3.6-anaconda-4.4

# Copy jupyter notebook in temp file
jupyter nbconvert --to notebook $inputipynb --output temp
# Comment instances of all jupyter magics
sed -i'' 's/%%time/#%%time/g' temp.ipynb
sed -i'' 's/%load_ext/#%load_ext/g' temp.ipynb
sed -i'' 's/%matplotlib/#%matplotlib/g' temp.ipynb
sed -i'' 's/%autoreload/#%autoreload/g' temp.ipynb
# Comment figure display
sed -i'' 's/plt.show()/plt.close()/g' temp.ipynb
# Uncomment matplotlib option to allow plotting on a non interactive backend
sed -i'' 's/# matplotlib.use(/matplotlib.use(/g' temp.ipynb

# Convert temp copy to python script
jupyter nbconvert --to script temp.ipynb --output $outputname

# remove temp files
rm temp.ipynb
