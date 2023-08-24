#!/bin/bash

###
### Script to generate .txt file table of scans for preprocessing pipeline
###

# Set original (Steve) data directory
orig_dir=/rds/project/rds-83zLpPieDV8/q10014/nii

# Set output directory
out_dir=/rds/project/rds-cAQcxgoLHoA/Livia

# Set study ID prefix
study='B3526'

# Get list of subjects
subjects=$(ls ${orig_dir} | grep ^${study}*)

# Create csv where column 1 is subject ID and column 2 is scan date
for sub in ${subjects}
do 
	for ses in $(ls ${orig_dir}/$sub)
	do 
		echo "$sub,$ses"
	done
done > ${out_dir}/scan_list.txt
