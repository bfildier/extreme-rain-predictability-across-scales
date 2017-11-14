#!/bin/bash
#SBATCH --mail-user=bfildier@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH --job-name="analyzePointWiseScalingOmega500TsPs"
#SBATCH --output="../logs/analyzePointWiseScalingOmega500TsPs.%j.%N.out"
##SBATCH --partition=regular
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00
#SBATCH --constraint=haswell
##SBATCH --constraint=knl,quad,cache

module load python/3.6-anaconda-4.4

# Convert notebook to python script
./ipynb2py36.sh analyzePointWiseScalingOmega500TsPs.ipynb analyzePointWiseScalingOmega500TsPs

# Edit options in script
startdate=185005010000
enddate=185005020000

sed -i='' 's/^dates =.*/dates = "'$startdate'","'$enddate'"/g' analyzePointWiseScalingOmega500TsPs.py

sed -i='' 's/^daskarray =.*/daskarray = False/' analyzePointWiseScalingOmega500TsPs.py


#---- Run script----#
# Start timer
start=`date +%s`
todaysdate=`date +"%Y%m%d-%H%M"`

python -u analyzePointWiseScalingOmega500TsPs.py &> ../logs/analyzePointWiseScalingOmega500TsPs_${todaysdate}.log

# End timer and save time execution
end=`date +%s`
runtime=$((end-start))
python -c "import datetime as dt; t_el=dt.timedelta(seconds=$runtime);\
print('Script successfully terminated in %s.'%t_el)"
