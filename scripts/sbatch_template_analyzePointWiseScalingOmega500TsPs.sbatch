#!/bin/bash
#SBATCH --mail-user=bfildier@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH --job-name="analyzePointWiseScalingOmega500TsPs"
#SBATCH --output="../logs/analyzePointWiseScalingOmega500TsPs.%j.%N.out"
#SBATCH --account=m1517
#SBATCH --partition=regular
##SBATCH --partition=shared
##SBATCH --partition=debug
#SBATCH --qos=premium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
##SBATCH --cpus-per-task=30
##SBATCH --time=00:30:00
#SBATCH --time=14:00:00
#SBATCH --constraint=haswell
##SBATCH --constraint=knl,quad,cache

module load python/3.6-anaconda-4.4

todaysdate=`date +"%Y%m%d-%H%M"`
scriptname=analyzePointWiseScalingOmega500TsPs_${todaysdate}_${SLURM_JOB_ID}

# Convert notebook to python script
./ipynb2py36.sh analyzePointWiseScalingOmega500TsPs.ipynb $scriptname 

# Edit options in script
startdate=185005010100
enddate=185105010000
#enddate=185005020000
sed -i='' 's/^dates =.*/dates = "'$startdate'","'$enddate'"/g' ${scriptname}.py
daskarray=True
sed -i='' 's/^daskarray =.*/daskarray = '${daskarray}'/'  ${scriptname}.py
tracktime=True
sed -i='' 's/^tracktime =.*/tracktime = '${tracktime}'/' ${scriptname}.py
compset=FAMIPC5
sed -i='' 's/^compset =.*/compset = "'${compset}'"/' ${scriptname}.py
experiment=abrupt4xCO2
sed -i='' 's/^experiment =.*/experiment = "'${experiment}'"/' ${scriptnam}.py

#---- Run script----#

python -u ${scriptname}.py > ../logs/${scriptname}.log

