#!/bin/bash
 
# Request resources:
#SBATCH --job-name=julia_cm_ntasks10
#SBATCH --output=julia_cm_%j.out
#SBATCH --error=julia_cm_%j.err
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1   # Number of CPU core
#SBATCH --mem=100G          # memory required, up to 250G on standard nodes.
#SBATCH --time=1:00:0       # time limit for job (format:  days-hours:minutes:seconds)
#SBATCH --gres=tmp:1G       # temporary disk space required on the compute node ($TMPDIR),
#                        up to 400G
# Run in the 'shared' queue (job may share node with other jobs)
#SBATCH -p shared
 
# Commands to be run:
module load julia

export JULIA_NUM_THREADS=1
export JULIA_BLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MALLOC_TRIM_THRESHOLD_=131072

echo "SLURM job $SLURM_JOB_ID on $(hostname)"
echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"

julia --project=. L1_Wolfe.jl 2 12 2 12 100