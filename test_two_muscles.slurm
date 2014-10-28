#!/bin/bash

#SBATCH --job-name=tytell2musc
#SBATCH --output=slurm.%N.%j.out
#SBATCH --error=slurm.%N.%j.err
#SBATCH --nodes=1

#SBATCH --partition=batch

#SBATCH --mail-user=eric.tytell@tufts.edu

# Finally, Start the program
module load matlab
matlab -nojvm -nodisplay < run_test_two_muscles.m > test_two_muscles.out