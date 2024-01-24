#!/bin/bash
#SBATCH -N 2
#SBATCH -n 96
#SBATCH -t 01:00:00
#SBATCH -p workq
#SBATCH -A loni_ceds3d
#SBATCH -o o.out
#SBATCH -e e.err
#SBATCH -J ship_in_current
#load proteus module and ensure proteus's python is in path

date

module purge
module load proteus/1.8.1

mkdir -p $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
cd "$SLURM_SUBMIT_DIR"
cp *.pyx $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cp *.py $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cp *.sh $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cp *.stl $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cd $WORK/$SLURM_JOB_NAME.$SLURM_JOBID

python setup.py build_ext -i
srun parun --TwoPhaseFlow ship_in_current_0.py -F -l 5 -C "he=0.08"

exit 0
