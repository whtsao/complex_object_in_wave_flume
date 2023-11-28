#!/bin/bash
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 00:10:00
#SBATCH -p workq
#SBATCH -A hpc_ceds3d
#SBATCH -o o.out
#SBATCH -e e.err
#SBATCH -J ship_in_current
#load proteus module and ensure proteus's python is in path

date

module purge
module load proteus/fct
module load intel/2021.5.0
module load mvapich2/2.3.7/intel-2021.5.0
module load gcc/11.2.0
export LD_LIBRARY_PATH=/home/packages/compilers/intel/compiler/2022.0.2/linux/compiler/lib/intel64_lin:${LD_LIBRARY_PATH}
export MV2_HOMOGENEOUS_CLUSTER=1

mkdir -p $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
cd $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
cp $SLURM_SUBMIT_DIR/petsc.options.superlu_dist .
cp $SLURM_SUBMIT_DIR/*.py .
cp $SLURM_SUBMIT_DIR/*.sh .

srun parun --TwoPhaseFlow ship.py -F -l 5 -C "he=0.1"

exit 0
