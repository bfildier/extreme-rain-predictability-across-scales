#!/bin/bash

./analyzeScalingAllScales.sh | tee logs/analyzeScaling_`date +"%Y%m%d-%H%M"`.log
