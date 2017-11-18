#!/bin/bash

module load python/3.6-anaconda-4.4

##-- Set working directory --##
SCRIPTDIR=automated_scripts
if [ ! -d "$SCRIPTDIR" ]
then
    mkdir $SCRIPTDIR
    mkdir $SCRIPTDIR/logs
fi 

##-- Run options--##
todaysdate=`date +"%Y%m%d-%H%M"`
template_nameroot=analyzeScaling_template_$todaysdate
template_analysis_script=${template_nameroot}.py
template_batch_script=sbatch_${template_nameroot}.sbatch

##-- Analysis script options --##
startdate=185005010100
enddate=185105010000
daskarray=False
tracktime=True
compsets='FSPCAMm_AMIP'
experiments='piControl'
time_strides='1hr'
resolutions='3dx'


##-- Batch script options --##
# WRITE HERE SBATCH OPTIONS

##-- Main --##

# Template scripts
./ipynb2py36.sh analyzePointWiseScalingOmega500TsPs_template.ipynb $SCRIPTDIR/${template_nameroot}
cp sbatch_template_analyzePointWiseScalingOmega500TsPs.sbatch $SCRIPTDIR/${template_batch_script}

cd $SCRIPTDIR

# Correct paths and change global options in analysis script
sed -i'' 's/(workdir)/(os.path.dirname(workdir))/g' ${template_analysis_script}
sed -i'' 's/^daskarray =.*/daskarray = '${daskarray}'/' ${template_analysis_script}
sed -i'' 's/^tracktime =.*/tracktime = '${tracktime}'/' ${template_analysis_script}


# Launch all runs
for compset in `echo ${compsets}`;
do
    for experiment in `echo ${experiments}`;
    do
        for time_stride in `echo ${time_strides}`;
        do
            for resolution in `echo ${resolutions}`;
            do

##-- Create scripts --##
                 
                nameroot=${template_nameroot}_${compset}_${experiment}_${time_stride}_${resolution}
                analysisscript=${nameroot}.py
                batchscript=sbatch_${nameroot}.sbatch

                # Duplicate analysis script
                cp ${template_analysis_script} ${analysisscript}
                # Duplicate batch script
                cp ${template_batch_script} $batchscript

##-- Analysis script options --#

                sed -i'' 's/^dates =.*/dates = "'$startdate'","'$enddate'"/' ${analysisscript}
                sed -i'' 's/^compset =.*/compset = "'${compset}'"/' ${analysisscript}
                sed -i'' 's/^experiment =.*/experiment = "'${experiment}'"/' ${analysisscript}
                sed -i'' 's/^time_stride =.*/time_stride = "'${time_stride}'"/' ${analysisscript}
                sed -i'' 's/^resolution =.*/resolution = "'${resolution}'"/' ${analysisscript}
                
##-- Batch script options --#
                
                sed -i'' 's/^#SBATCH --output=.*/#SBATCH --output="logs\/'${nameroot}'.%j.%N.out"/' $batchscript
                sed -i'' 's/^scriptname=.*/scriptname='${nameroot}'/' $batchscript

##-- Launch batch script--##

                sbatch $batchscript
                
            done
        done
    done
done


exit 0
