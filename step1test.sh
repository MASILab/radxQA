#!/bin/bash
for i in $(find /blazepool/radx/QA/NEWDATA/* -type d -mindepth 1  -maxdepth 1); do
	# ** may have issues with spaces in file names ** 
	# for i in /GRAID/BlazeFake/NEWDATA/
	# echo $i;
	# only if EPI_STABILITY or T1_axial or ACR_T1
	
	if [[ $i == *"EPI_STABILITY" ]] || [[ $i == *"T1_axial" ]] || [[ $i == *"ACR" ]]; 
	then
		echo $i;
		
	fi
done