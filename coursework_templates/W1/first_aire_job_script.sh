#!/bin/bash
#SBATCH --job-name=first_aire_job
#SBATCH --output=output_%j.out        # Output file (%j = job ID)
#SBATCH --error=error_%j.err          # Error file (%j = job ID)
#SBATCH --time=00:01:00
#SBATCH --mem=1G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# Run the job
python ~/COMP3270_work/first_aire_job.py

