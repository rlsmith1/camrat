#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-825
#SBATCH --cpus-per-task=1
#SBATCH -e /rds/project/rds-cAQcxgoLHoA/afni_preprocessing/scripts/pipeline/logs/03.roi-stats.err
#SBATCH -o /rds/project/rds-cAQcxgoLHoA/afni_preprocessing/scripts/pipeline/logs/03.roi-stats.out
#SBATCH -t 00:05:00
#SBATCH --mail-type=ALL

export PATH=$PATH:/rds/project/rds-cAQcxgoLHoA/abin-lena/

base_dir=/rds/project/rds-cAQcxgoLHoA/afni_preprocessing
data_dir=${base_dir}/data/derived
outdir=${base_dir}/outputs/stats

scan_list=($(cat ${base_dir}/scan_list.txt))
export line=${scan_list[${SLURM_ARRAY_TASK_ID}]}

IFS=',' read sub ses cont init_scale <<< ${line}

echo ${sub}
echo ${ses}
echo ${cont}

scan_dir=${data_dir}/${sub}/${ses}/anat/02.warp_struct_${cont}
infile=${scan_dir}/${sub}_${ses}_${cont}_deobliqued_RIP_shft.nii.gz
mset=${scan_dir}/SIGMA_Anatomical_Brain_Atlas_in_${sub}_${ses}_${cont}_deobliqued_RIP_shft.nii.gz

if [[ -f ${infile} ]]; then

        echo "Calculating ROI stats for "${infile}
	3dROIstats -mask ${mset} -nzvolume ${infile} > ${outdir}/${sub}_${ses}_${cont}_3dROIstats.txt

else
        echo "Can't find "${infile}
fi

