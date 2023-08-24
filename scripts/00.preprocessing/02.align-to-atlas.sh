#!/bin/bash
#SBATCH -p cclake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -e /rds/project/rds-cAQcxgoLHoA/Livia/scripts/logs/02.err
#SBATCH -o /rds/project/rds-cAQcxgoLHoA/Livia/scripts/logs/02.out
#SBATCH -t 00:05:00
#SBATCH --mail-type=ALL

# point path to updated AFNI
export PATH=$PATH:/rds/project/rds-cAQcxgoLHoA/abin-lena/

# set paths
base_dir=/rds/project/rds-cAQcxgoLHoA/Livia
data_dir=${base_dir}/data
derived_dir=${base_dir}/data/derived
template=${base_dir}/SIGMA_v1.2/SIGMA_Rat_Anatomical_Imaging/SIGMA_Rat_Anatomical_ExVivo_Template/SIGMA_ExVivo_Brain_Template.nii

mkdir -p ${derived_dir}

#subject [sub], session [ses] and contrast [cont] selected in parent shell script (run-02.sh)
sub=$1
ses=$2
cont=$3

# set infile (native scan)
infile=${data_dir}/${sub}/${ses}/anat/${sub}_${ses}_${cont}.nii

# create corresponding directories in derived_dir
out_dir=${derived_dir}/${sub}/${ses}/anat
mkdir -p ${out_dir}

# deoblique
if [[ -f ${infile} ]]; then
	
	echo "preprocessing scan "${infile}
	
        is_oblique=$( 3dinfo -is_oblique ${infile} )
        if [[ ${is_oblique} ]]; then
		3dWarp -oblique2card -prefix ${out_dir}/${sub}_${ses}_${cont}_deobliqued.nii ${infile}
        else
             	cp ${infile} ${out_dir}/${sub}_${ses}_${cont}_deobliqued.nii
        fi
fi


# rotate scan to align orientation with SIGMA template
3dresample -orient RAI -prefix ${out_dir}/${sub}_${ses}_${cont}_deobliqued_RAI.nii -input ${out_dir}/${sub}_${ses}_${cont}_deobliqued.nii
3dcopy ${out_dir}/${sub}_${ses}_${cont}_deobliqued_RAI.nii ${out_dir}/${sub}_${ses}_${cont}_deobliqued_RIP.nii
3drefit -orient RIP ${out_dir}/${sub}_${ses}_${cont}_deobliqued_RIP.nii

# Center of mass align
@Align_Centers -base ${template} -dset ${out_dir}/${sub}_${ses}_${cont}_deobliqued_RIP.nii
