#!/bin/bash

./analyzeScalingAllScalesMerged.sh | tee logs/analyzeScalingMerged_`date +"%Y%m%d-%H%M"`.log
