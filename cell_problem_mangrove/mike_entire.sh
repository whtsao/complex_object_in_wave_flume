#!/bin/bash
#SBATCH -N 4
#SBATCH -n 256
#SBATCH -t 01:00:00
#SBATCH -p workq
#SBATCH -A hpc_ceds3d
#SBATCH -J wave_pass_realtree
#SBATCH -o o.out
#SBATCH -e e.err
#load proteus module and ensure proteus's python is in path
date
module load proteus/fct
export LD_LIBRARY_PATH=/home/packages/compilers/intel/compiler/2022.0.2/linux/compiler/lib/intel64_lin:${LD_LIBRARY_PATH}
export MV2_HOMOGENEOUS_CLUSTER=1
mkdir -p $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
cd $SLURM_SUBMIT_DIR
cp setup.py $WORK/$SLURM_JOB_NAME.$SLURM_JOBID 
cp cell.py $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cp mangrove.pyx $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cp mangrove_2.pyx $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cp nse_p.py $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cp nse_n.py $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cp petsc.options.asm $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cp *.sh $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
cd $WORK/$SLURM_JOB_NAME.$SLURM_JOBID
python setup.py build_ext -i
#he = .05 for 1 core, .025 8 core, .0125 64 core
srun parun nse_p.py nse_n.py -F -l 5 -C "he=.02 usePETSc=True cell_bottom=0. cell_height=1. T=1." -O petsc.options.asm -D p140
exit 0
