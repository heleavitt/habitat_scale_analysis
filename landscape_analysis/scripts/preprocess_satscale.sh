#!/bin/bash
#SBATCH --job-name=pre_processing
#SBATCH --output=pre_processing.out
#SBATCH --error=pre_processing.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=3

cd $SLURM_SUBMIT_DIR

ml Python/3.9.6-GCCcore-11.2.0
. ~/env/fourchonenv/bin/activate

python pre_processing.py

