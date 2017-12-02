#!/bin/bash

#experiments_all="'piControl' 'abrupt4xCO2'"
experiments_all="'piControl'"
#compsets_all="'FSPCAMm_AMIP' 'FAMIPC5'"
compsets_all="'FAMIPC5'"
time_strides_all="'1h3h6h12h' '1d2d4d8d'"
script='analyzeScalingAllScales.sh'

for experiments in `echo ${experiments_all}`; do
    
    for compsets in `echo ${compsets_all}`; do
    
        for time_strides in `echo ${time_strides_all}`; do
        
            ts=`echo ${time_strides} | sed 's/h/h\ /g' | sed 's/d/d\ /g'`
            
            date
            echo
            echo ${experiments} ${compsets} ${ts}
            
            sed -i'' 's/^experiments=.*/experiments='${experiments}'/' $script
            sed -i'' 's/^compsets=.*/compsets='${compsets}'/' $script
            sed -i'' "s/^time_strides=.*/time_strides=${ts}/" $script

            ./launch_analysis.sh
            echo
            echo "wait 5 hours before next batch submissions"
            echo

            sleep 18000
        
        done
    done
done

exit 0
