#!/bin/bash
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 24:00:00
#SBATCH -p workq
#SBATCH -A hpc_ceds3d
#SBATCH -J cell
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
#srun parun nse_p.py nse_n.py -F -l 5 -C "he=.0125 usePETSc=True cell_bottom=0.0 cell_height=0.05 T=100.0" -O petsc.options.asm  -D p00
#srun parun nse_p.py nse_n.py -F -l 5 -C "he=.0125 usePETSc=True cell_bottom=0.05 cell_height=0.05 T=100.0" -O petsc.options.asm -D p05
#srun parun nse_p.py nse_n.py -F -l 5 -C "he=.0125 usePETSc=True cell_bottom=0.10 cell_height=0.05 T=100.0" -O petsc.options.asm -D p10
#srun parun nse_p.py nse_n.py -F -l 5 -C "he=.0125 usePETSc=True cell_bottom=0.15 cell_height=0.05 T=100.0" -O petsc.options.asm -D p15
#srun parun nse_p.py nse_n.py -F -l 5 -C "he=.0125 usePETSc=True cell_bottom=0.20 cell_height=0.05 T=100.0" -O petsc.options.asm -D p20
#srun parun nse_p.py nse_n.py -F -l 5 -C "he=.0125 usePETSc=True cell_bottom=0.25 cell_height=0.05 T=100.0" -O petsc.options.asm -D p25
#srun parun nse_p.py nse_n.py -F -l 5 -C "he=.0125 usePETSc=True cell_bottom=0.30 cell_height=0.05 T=100.0" -O petsc.options.asm -D p30
#srun parun nse_p.py nse_n.py -F -l 5 -C "he=.0125 usePETSc=True cell_bottom=0.35 cell_height=0.05 T=100.0" -O petsc.options.asm -D p35
#srun parun nse_p.py nse_n.py -F -l 5 -C "he=.0125 usePETSc=True cell_bottom=0.40 cell_height=0.05 T=100.0" -O petsc.options.asm -D p40
#srun parun nse_p.py nse_n.py -F -l 5 -C "he=.0125 usePETSc=True cell_bottom=0.45 cell_height=0.05 T=100.0" -O petsc.options.asm -D p45
srun parun nse_p.py nse_n.py -F -l 5 -C "he=.0125 usePETSc=True cell_bottom=0.50 cell_height=0.05 T=100.0" -O petsc.options.asm -D p50
exit 0
