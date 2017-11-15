#!/bin/bash
#SBATCH --mail-user=bfildier@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH --job-name="analyzePointWiseScalingOmega500TsPs"
#SBATCH --output="../logs/analyzePointWiseScalingOmega500TsPs.%j.%N.out"
#SBATCH --account=m1517
##SBATCH --partition=regular
##SBATCH --partition=shared
#SBATCH --partition=debug
##SBATCH --qos=premium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
##SBATCH --cpus-per-task=30
#SBATCH --time=00:30:00
##SBATCH --time=9:00:00
#SBATCH --constraint=haswell
##SBATCH --constraint=knl,quad,cache

module load python/3.6-anaconda-4.4

todaysdate=`date +"%Y%m%d-%H%M"`
scriptname=analyzePointWiseScalingOmega500TsPs_${todaysdate}_${SLURM_JOB_ID}

# Convert notebook to python script
./ipynb2py36.sh analyzePointWiseScalingOmega500TsPs.ipynb $scriptname 

# Edit options in script
startdate=185005010100
enddate=185005020000
#enddate=185005020000
sed -i='' 's/^dates =.*/dates = "'$startdate'","'$enddate'"/g' ${scriptname}.py
daskarray=False
sed -i='' 's/^daskarray =.*/daskarray = '${daskarray}'/'  ${scriptname}.py
tracktime=True
sed -i='' 's/^tracktime =.*/tracktime = '${tracktime}'/' ${scriptname}.py
compset=FSPCAMm_AMIP
sed -i='' 's/^compset =.*/compset = "'${compset}'"/' ${scriptname}.py
experiment=piControl
sed -i='' 's/^experiment =.*/experiment = "'${experiment}'"/' ${scriptname}.py
time_stride='1hr'
sed -i='' 's/^time_stride =.*/time_stride = "'${time_stride}'"/' ${scriptname}.py

#---- Run script----#

python -u ${scriptname}.py > ../logs/${scriptname}.log

