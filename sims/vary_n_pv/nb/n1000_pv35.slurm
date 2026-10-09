#!/bin/bash
#SBATCH --mail-user=lotspes@wfu.edu  # email address
#SBATCH --mail-type=BEGIN  # Alerts sent when job begins, ends, or aborts
#SBATCH --nodes=1   # Number of nodes required
#SBATCH --ntasks=1   # Number of nodes required
#SBATCH --mem=116G  # Total memory (RAM) required, per node
#SBATCH --time=0-22:30:00  # Wall Clock time (dd-hh:mm:ss) [max of 14 days]
#SBATCH --array=0-19
#SBATCH --job-name=1000/35 # job name
#SBATCH --partition=small

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load apps/r/4.3.3 
Rscript n1000_pv35.R $SLURM_ARRAY_TASK_ID
