#!/bin/bash
#SBATCH -p cclake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -e /rds/project/rds-cAQcxgoLHoA/Livia/scripts/logs/03.err
#SBATCH -o /rds/project/rds-cAQcxgoLHoA/Livia/scripts/logs/03.out
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL

# point path to updated AFNI
export PATH=$PATH:/rds/project/rds-cAQcxgoLHoA/abin-lena/

# set paths to atlas, template, and data
base_dir=/rds/project/rds-cAQcxgoLHoA/Livia
template=${base_dir}/SIGMA_v1.2/SIGMA_Rat_Anatomical_Imaging/SIGMA_Rat_Anatomical_ExVivo_Template/SIGMA_ExVivo_Brain_Template.nii
atlas=${base_dir}/SIGMA_v1.2/SIGMA_Rat_Brain_Atlases/SIGMA_Anatomical_Atlas/SIGMA_Anatomical_Brain_Atlas.nii
mask=${base_dir}/SIGMA_v1.2/SIGMA_Rat_Anatomical_Imaging/SIGMA_Rat_Anatomical_ExVivo_Template/SIGMA_ExVivo_Brain_Mask.nii

data_dir=${base_dir}/data/derived

# session, contrast, and subject are defined in the 'parent' shell script (run-03.sh)
sub=$1
ses=$2
cont=$3
init_scale=$4

# set directory based on these parameters
subj_dir=${data_dir}/${sub}/${ses}/anat
infile=${subj_dir}/${sub}_${ses}_${cont}_deobliqued_RIP_shft.nii

if [[ -f ${infile} ]]; then
        echo "Warping file "${infile}
        @animal_warper \
                -input ${infile} \
                -base  ${template} \
                -atlas  ${atlas} \
                -skullstrip ${mask} \
                -outdir ${subj_dir}/03.register_${cont} \
		-cost lpa+ZZ \
		-init_scale ${init_scale}
else
        echo "Can't find "${infile}
fi
