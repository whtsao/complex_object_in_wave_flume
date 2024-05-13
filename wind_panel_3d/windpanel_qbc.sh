#!/bin/bash
#SBATCH -N 32
#SBATCH -n 1536
#SBATCH -t 24:00:00
#SBATCH -p workq
#SBATCH -A loni_ceds3d
#SBATCH -J panel_in_wind
#SBATCH -o o.out
#SBATCH -e e.err
#load proteus module and ensure proteus's python is in path

date
module purge
module load proteus/1.8.1

cd $SLURM_SUBMIT_DIR
cp *.py $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
cp *.pyx $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cp *.asm $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cp *.sh $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cp *.stl $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cp *.csv $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cd $WORK/$SLURM_JOB_NAME.$SLURM_JOBID

python setup.py build_ext -i
#he = .05 for 1 core, .025 8 core, .0125 64 core

# use -N 16 -n 1024 -t 72 for this case
#srun parun nse_p.py nse_n.py -F -l 5 -C "he=0.04" -O petsc.options.asm #-D pcube
srun parun nse_p.py nse_n.py -F -l 5 -C "he=0.06" -O petsc.options.asm #-D pcube

#srun parun nse_p.py nse_n.py -F -l 5 -C "he=0.1" -O petsc.options.asm #-D pcube

exit 0

