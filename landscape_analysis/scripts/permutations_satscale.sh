#!/bin/bash
#SBATCH --job-name=process_perm
#SBATCH --output=logs/process_perm_%A_%a.out
#SBATCH --error=logs/process_perm_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --array=0-35

cd $SLURM_SUBMIT_DIR

ml Python/3.9.6-GCCcore-11.2.0
. ~/env/fourchonenv/bin/activate

# Define edge_distances and buffer_distances arrays
edge_distances=(1 3 5)
buffer_distances=(100 150 200 250 300 400 500 600 700 800 900 1000)

# Calculate total number of permutations
num_edge=${#edge_distances[@]}
num_buffer=${#buffer_distances[@]}
total_permutations=$((num_edge * num_buffer))

# Get the SLURM array task ID
task_id=$SLURM_ARRAY_TASK_ID

# Calculate the indices for edge_distance and buffer_distance
edge_index=$(( task_id / num_buffer ))
buffer_index=$(( task_id % num_buffer ))

# Get the edge_distance and buffer_distance
edge_distance=${edge_distances[$edge_index]}
buffer_distance=${buffer_distances[$buffer_index]}

echo "Processing edge_distance=${edge_distance}, buffer_distance=${buffer_distance}"

python scripts/satscale_permutations_241022.py $edge_distance $buffer_distance
