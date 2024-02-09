#!/bin/bash
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -t 00:10:00
#SBATCH -p workq
#SBATCH -A loni_ceds3d
#SBATCH -o o.out
#SBATCH -e e.err
#SBATCH -J floating_panels_2d
#load proteus module and ensure proteus's python is in path

date

module purge
module load proteus/1.8.1

mkdir -p $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
cd $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
cp $SLURM_SUBMIT_DIR/petsc.options.superlu_dist .
cp $SLURM_SUBMIT_DIR/*.py .
cp $SLURM_SUBMIT_DIR/*.sh .

srun parun --TwoPhaseFlow floating_panels.py -F -l 5 -C "he=0.2" # all on
#srun parun --TwoPhaseFlow floating_panels.py -F -l 5 -C "he=0.2 bodybool2=False linkbool=False" # only one panel
#srun parun --TwoPhaseFlow floating_panels.py -F -l 5 -C "he=0.2 linkbool=False" # two panels unconnected

exit 0
