#!/bin/bash
 
# Request resources:
#SBATCH --job-name=LinfFM
#SBATCH --output=L_{inf}_13_13_1000/LinfFM_%j.out
#SBATCH --error=L_{inf}_13_13_1000/LinfFM_%j.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1   # Number of CPU core
#SBATCH --mem=20G           # memory required, up to 250G on standard nodes.
#SBATCH --time=0:45:0       # time limit for job (format:  days-hours:minutes:seconds)
#SBATCH --gres=tmp:1G       # temporary disk space required on the compute node ($TMPDIR),
#                        up to 400G
# Run in the 'shared' queue (job may share node with other jobs)
#SBATCH -p shared

#SBATCH --array=1-1000
 
# Commands to be run:
module load julia

export JULIA_NUM_THREADS=1
export JULIA_BLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

echo "SLURM job $SLURM_JOB_ID on $(hostname)"
echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"

julia --project=. Linf_Wolfe.jl 2 13 2 13 ${SLURM_ARRAY_TASK_ID} ${SLURM_ARRAY_TASK_ID}