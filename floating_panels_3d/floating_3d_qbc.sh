#!/bin/bash
#SBATCH -N 8
#SBATCH -n 384
#SBATCH -t 00:20:00
#SBATCH -p workq
#SBATCH -A loni_ceds3d
#SBATCH -o o.out
#SBATCH -e e.err
#SBATCH -J floating_panels_3d
#load proteus module and ensure proteus's python is in path

date

module purge
module load proteus/1.8.1

mkdir -p $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
cd $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
cp $SLURM_SUBMIT_DIR/petsc.options.superlu_dist .
cp $SLURM_SUBMIT_DIR/*.py .
cp $SLURM_SUBMIT_DIR/*.sh .
cp $SLURM_SUBMIT_DIR/*.stl .

srun parun TN_with_box_so.py -F -l 5 -C "he=2. T=5."

# 1: one panel
#srun parun --TwoPhaseFlow floating_panels.py -F -l 5 -C "he=0.1 bodybool2=False linkbool=False mooring=False"

# 2: two panels
#srun parun --TwoPhaseFlow floating_panels.py -F -l 5 -C "he=0.1 bodybool2=True linkbool=False mooring=False"

# 3: one panel moored
#srun parun --TwoPhaseFlow floating_panels.py -F -l 5 -C "he=0.1 bodybool2=False linkbool=False mooring=True"

# 4: two panel moored without link
#srun parun --TwoPhaseFlow floating_panels.py -F -l 5 -C "he=0.1 bodybool2=True linkbool=False mooring=True"

# 5: two panel linked without mooring
#srun parun --TwoPhaseFlow floating_panels.py -F -l 5 -C "he=0.1 bodybool2=True linkbool=True mooring=False"

# 6: two panel linked without mooring
#srun parun --TwoPhaseFlow floating_panels.py -F -l 5 -C "he=0.1 bodybool2=True linkbool=True mooring=True"

exit 0
