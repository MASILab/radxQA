#!/bin/bash
# loop through all folders in NEWDATA directory
# run run_dcmdjpeg.sh on each folder in NEWDATA
# run mcverter on each 
# move nii and txt to NIFTIS
# delete folder in NEWDATA

# loop through all folders in NEWDATA directory
# for i in $(find $PWD/NEWDATA -maxdepth 1 -type d); 
for i in $(find /Volumes/GRAID/BlazeFake/QA/NEWDATA/* -type d -mindepth 1  -maxdepth 1); do
	# ** may have issues with spaces in file names ** 
	# for i in /GRAID/BlazeFake/NEWDATA/
	echo $i;
	# only if EPI_STABILITY or T1_axial or ACR_T1
	if [[ $i == *"EPI_STABILITY" ]] || [[ $i == *"T1_axial" ]] || [[ $i == *"ACR" ]]; then
	
		if test -d $i/decompressed_dicoms/OUT; 
		then echo 'Already Done...'; 
		else
			# run run_dcmdjpeg.sh on each folder in NEWDATA 
			echo "running dcmdjpeg..."
			/Volumes/GRAID/BlazeFake/QA/EXTRA/run_dcmdjpeg.sh $i
			# run mcverter on each
			echo "running mcverter..."
			mcverter -o $i/decompressed_dicoms/OUT/ -d -f nifti -n -F +PatientName+PatientId+SeriesDate+SeriesTime+StudyId+StudyDescription+SeriesNumber+SequenceName+ProtocolName+SeriesDescription -v $i/decompressed_dicoms/
			# move  nii and txt to NIFTIs
			echo "moving data to NIFTIS..."
			mv $i/decompressed_dicoms/OUT/* /Volumes/GRAID/BlazeFake/QA/NIFTIS
			# delete folder in NEWEDATA
			echo "deleting folder in NEWDATA..."
			# rm -rf $i
		fi
	
	
	
	fi
done


