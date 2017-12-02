#!/bin/bash

scriptin=$1
scriptout=$2

sed -e '/# start comment/,/# stop comment/ s/^/# /' $scriptin > temp1
sed -e '/# # start uncomment/,/# # stop uncomment/ s/^# //' temp1 > temp2
sed -e '/#++#/,/#--#/ s/^/        /' temp2 > temp3

mv temp3 $scriptout
rm temp*