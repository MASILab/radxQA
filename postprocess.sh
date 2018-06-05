#!/bin/bash
# loop through all folders in PROCESSED directory
# find all .sh files
# run them all (calls to python phantom2redcap.py)
# delete .sh files

for i in $(find /Volumes/GRAID/BlazeFake/QA/PROCESSED/*.sh -type f); do
	echo $i;
	# run them all
	echo "running SH script..."
	sh $i
	# delete them
	echo "deleting SH script..."
	rm -rf $i
done

