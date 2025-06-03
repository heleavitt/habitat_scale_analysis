#!/bin/bash
#SBATCH --job-name=drone_landscap
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --cpus-per-task=50
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=280gb                     # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=logs/%x.%j.out            # Standard output log
#SBATCH --error=logs/%x.%j.err             # Standard error log

#SBATCH --mail-type=ALL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=hl51981@uga.edu  # Where to send mail	

cd $SLURM_SUBMIT_DIR
ml Python/3.9.6-GCCcore-11.2.0
. ~/env/fourchonenv/bin/activate

python scripts/smallscale_permutations_241029.py
