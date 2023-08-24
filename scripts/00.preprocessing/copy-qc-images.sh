#!/bin/bash

# Set paths to directories where data are
base_dir=/rds/project/rds-cAQcxgoLHoA/Livia
derived_dir=${base_dir}/data/derived

# Set output directory (where you want QC images to go
out_dir=${base_dir}/outputs/QC

# Read scan table for list of scans
scan_table=$(cat ${base_dir}/scan_table_for_registration.txt)

# Loop through scans
for row in ${scan_table}
do

	# Pull subject, session, and contrast info from table
	IFS=',' read sub ses cont init_scale <<< ${row}

	echo ${sub}
	echo ${ses}
	echo ${cont}

	# Locate QC directory
	qc_dir=${derived_dir}/${sub}/${ses}/anat/03.register_${cont}/QC

	# Identify QC images of interest (in this case, qc_03)
	images=${qc_dir}/qc_03.input+wrpd_ATL_SIGMA_Anatomical_Brain_Atlas.*.png

	# Rename QC images with sub, ses, and cont, and copy to output directory
	for i in ${images}
	do 
		plane=$(echo $(basename $i) | cut -d '.' -f 3)
		cp $i ${out_dir}/${sub}_${ses}_${cont}_${plane}.png
	done

done

echo "finished"
