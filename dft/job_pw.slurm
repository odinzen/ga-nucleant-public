#!/bin/bash
#SBATCH --job-name=pw
#SBATCH --output=slurm-%x-%j.out
#SBATCH --partition=REPLACE_PARTITION      # e.g. general / htc -- set for Sol
#SBATCH --account=REPLACE_ACCOUNT
#SBATCH --nodes=1
#SBATCH --ntasks=32                        # match a node; tune per cell size
#SBATCH --time=08:00:00
#SBATCH --mem=0

# --- environment: load the Quantum ESPRESSO your site provides ---
module purge
module load quantum-espresso            # REPLACE with Sol's exact QE module name
# (or: module load quantum_espresso/7.x ; or source your own QE build)

# pw.x reads ./pw.in in the job's directory and writes pw.out
srun pw.x -in pw.in > pw.out

# tip: bulk_* are cheap (run first to get relaxed lattices); interface is the heavy one.
