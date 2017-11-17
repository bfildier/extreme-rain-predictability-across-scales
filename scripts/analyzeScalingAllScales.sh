#!/bin/bash

module load python/3.6-anaconda-4.4

## Directory to store all scripts
SCRIPTDIR=automated_scripts
if [ ! -d "$SCRIPTDIR" ]
then
    mkdir $SCRIPTDIR
    mkdir $SCRIPTDIR/logs
fi 

## Run options
todaysdate=`date +"%Y%m%d-%H%M"`

# Edit options in script
startdate=185005010100
enddate=185005020000
#enddate=185005020000
daskarray=False
tracktime=True
compsets='FSPCAMm_AMIP'
experiments='piControl'
time_strides='1hr'
resolutions='dx'


for compset in `echo ${compsets}`;
do
    for experiment in `echo ${experiments}`;
    do
        for time_stride in `echo ${time_strides}`;
        do
            for resolution in `echo ${resolutions}`;
            do

##-- Create scripts --##
                 
                nameroot=analyzePointWiseScalingOmega500TsPs_${todaysdate}_${compset}_${experiment}_${time_stride}_${resolution}
                analysisscript=${SCRIPTDIR}/${nameroot}.py
                batchscript=${SCRIPTDIR}/sbatch_${nameroot}.sbatch

                # Convert notebook to python script
                ./ipynb2py36.sh analyzePointWiseScalingOmega500TsPs_template.ipynb ${SCRIPTDIR}/${nameroot}
                # Duplicate batch script
                cp sbatch_template_analyzePointWiseScalingOmega500TsPs.sbatch $batchscript

##-- Python script options --#

                sed -i='' 's/^dates =.*/dates = "'$startdate'","'$enddate'"/g' ${analysisscript}
                sed -i='' 's/^daskarray =.*/daskarray = '${daskarray}'/' ${analysisscript}
                sed -i='' 's/^tracktime =.*/tracktime = '${tracktime}'/' ${analysisscript}
                sed -i='' 's/^compset =.*/compset = "'${compset}'"/' ${analysisscript}
                sed -i='' 's/^experiment =.*/experiment = "'${experiment}'"/' ${analysisscript}
                sed -i='' 's/^time_stride =.*/time_stride = "'${time_stride}'"/' ${analysisscript}
                sed -i='' 's/^resolution =.*/resolution = "'${resolution}'"/' ${analysisscript}
                
##-- Batch script options --#
                
                sed -i='' 's/^#SBATCH --output=.*/#SBATCH --output="../logs/'${scriptname}'.%j.%N.out"/' $batchscript
                sed -i='' 's/^scriptname=.*/scriptname='${nameroot}'/' $batchscript

##-- Send batch script to slurm --##

                sbatch $batchscript
                echo $batchscript" submitted"
                
            done
        done
    done
done

exit 0
