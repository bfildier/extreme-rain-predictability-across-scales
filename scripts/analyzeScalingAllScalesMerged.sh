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
#enddate=185005020000
daskarray=False
tracktime=True
bootstrap=True
plotAll2dPDFs=True
compsets='FSPCAMm_AMIP'
#compsets='FAMIPC5'
experiments='piControl'
#experiments='piControl abrupt4xCO2'
time_strides_all="\'1d\',\'3h\' \'6h\',\'12h\' \'1d\',\'2d\' \'4d\',\'8d\'"
#time_strides="\'1h\',\'3h\',\'6h\',\'12h\',\'1d\',\'2d\',\'4d\',\'8d\'"
#resolutions="\'1dx\',"
resolutions="\'1dx\',\'2dx\',\'3dx\',\'4dx\',\'5dx\',\'6dx\',\'7dx\',\'8dx\',\'9dx\'"

##-- Batch script options --##
runmode="regular"
#runmode="debug"

##-- Main --##

# Template scripts
./ipynb2py36.sh analyzePointWiseScaling_template.ipynb $SCRIPTDIR/${template_nameroot}
./reformat_script_for_multiple_scales.sh $SCRIPTDIR/${template_analysis_script} $SCRIPTDIR/${template_analysis_script}

echo $SCRIPTDIR/${template_analysis_script}

# batch ecript
cp sbatch_analyzePointWiseScaling_template.sbatch $SCRIPTDIR/${template_batch_script}


cd $SCRIPTDIR

# Correct paths and change global options in analysis script
sed -i'' 's/(workdir)/(os.path.dirname(workdir))/g' ${template_analysis_script}
sed -i'' 's/^daskarray =.*/daskarray = '${daskarray}'/' ${template_analysis_script}
sed -i'' 's/^tracktime =.*/tracktime = '${tracktime}'/' ${template_analysis_script}
sed -i'' 's/^bootstrap =.*/bootstrap = '${bootstrap}'/' ${template_analysis_script}
sed -i'' 's/^plotAll2dPDFs =.*/plotAll2dPDFs = '${plotAll2dPDFs}'/' ${template_analysis_script}
sed -i'' "s/^resolutions =.*/resolutions = ${resolutions}/" ${template_analysis_script}


# Change run options in batch script
if [ "$runmode" == "debug" ];
then
    sed -i'' 's/^#SBATCH --partition=.*/#SBATCH --partition=debug/' ${template_batch_script}
    sed -i'' 's/^#SBATCH --qos=/##SBATCH --qos=/' ${template_batch_script}
    sed -i'' 's/^#SBATCH --time=.*/#SBATCH --time=00:30:00/' ${template_batch_script}
    sed -i'' 's/^#SBATCH --mail-type=.*/#SBATCH --mail-type=ALL/' ${template_batch_script}
elif [ "$runmode" == "regular" ];
then 
    sed -i'' 's/^#SBATCH --partition=.*/#SBATCH --partition=regular/' ${template_batch_script}
    sed -i'' 's/.*#SBATCH --qos=.*/#SBATCH --qos=premium/' ${template_batch_script}
    sed -i'' 's/^#SBATCH --time=.*/#SBATCH --time=4:30:00/' ${template_batch_script}
    sed -i'' 's/^#SBATCH --mail-type=.*/#SBATCH --mail-type=FAIL,TIME_LIMIT_90/' ${template_batch_script}
fi

# Launch all runs
for compset in `echo ${compsets}`;
do
    for experiment in `echo ${experiments}`;
    do
        for time_strides in `echo ${time_strides_all}`;
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
            sed -i'' "s/^time_strides =.*/time_strides = ${time_strides}/" ${analysisscript}
    
##-- Batch script options --#
                    
            sed -i'' 's/^#SBATCH --output=.*/#SBATCH --output="logs\/'${nameroot}'.%j.%N.out"/' $batchscript
            sed -i'' 's/^scriptname=.*/scriptname='${nameroot}'/' $batchscript
    
##-- Launch batch script--##
    
            echo -n "$compset $experiment ${time_stride} ${resolution}: "
            sbatch $batchscript
    
        done
    done
done


exit 0
