#!/bin/bash
#SBATCH -p cclake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-159
#SBATCH -e /rds/project/rds-cAQcxgoLHoA/Livia/scripts/logs/run-02.err
#SBATCH -o /rds/project/rds-cAQcxgoLHoA/Livia/scripts/logs/run-02.out
#SBATCH -t 00:20:00
#SBATCH --mail-type=ALL

base_dir=/rds/project/rds-cAQcxgoLHoA/Livia
script_dir=${base_dir}/scripts

scan_table=($(cat ${base_dir}/scan_table_for_registration.txt))
export row=${scan_table[${SLURM_ARRAY_TASK_ID}]}

IFS=',' read sub ses cont init_scale <<< ${row}

echo ${sub}
echo ${ses}
echo ${cont}

# Run alignment script using values from table
bash ${script_dir}/02.align-to-atlas.sh ${sub} ${ses} ${cont}
