#!/bin/bash

remote_port=$1
local_port=$2

ssh cori.nersc.gov <<EOF
for p in `ps aux | grep ${remote_port} | grep -v grep | tr -s " " | cut -d" " -f2`; do
	echo $p
	kill $p
done
EOF

for p in `ps aux | grep ${local_port} | grep -v grep | tr -s " " | cut -d" " -f2`;
do
	echo $p
	kill $p
done

exit 0