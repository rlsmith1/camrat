#!/bin/bash
#SBATCH -p cclake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --array=0-159
#SBATCH -e /rds/project/rds-cAQcxgoLHoA/Livia/scripts/logs/run-03.err
#SBATCH -o /rds/project/rds-cAQcxgoLHoA/Livia/scripts/logs/run-03.out
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL

base_dir=/rds/project/rds-cAQcxgoLHoA/Livia
script_dir=${base_dir}/scripts

scan_table=($(cat ${base_dir}/scan_table_for_registration.txt))
export row=${scan_table[${SLURM_ARRAY_TASK_ID}]}

IFS=',' read sub ses cont init_scale <<< ${row}

echo ${sub}
echo ${ses}
echo ${cont}
echo ${init_scale}

# assign directory for @animal_warper outputs
out_dir=${base_dir}/data/derived/${sub}/${ses}/anat/03.register_${cont}

# erase what we had before, if it exists
if [[ -d ${out_dir} ]]; then
	rm -r ${out_dir}
fi

# run animal warper on scans using specified settings
bash ${script_dir}/03.register-to-atlas.sh ${sub} ${ses} ${cont} ${init_scale}
