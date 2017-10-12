#!/bin/bash

ssh cori.nersc.gov <<EOF
for p in `ps aux | grep 8889 | grep -v grep | tr -s " " | cut -d" " -f2`; do
	echo $p
	kill $p
done
EOF

for p in `ps aux | grep 8889 | grep -v grep | tr -s " " | cut -d" " -f2`; do
	echo $p
	kill $p
done