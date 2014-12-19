#!/bin/bash

#SBATCH --job-name=tytell2musc
#SBATCH --output=test_two_muscles.slurm.%N.%j.out
#SBATCH --error=test_two_muscles.slurm.%N.%j.err
#SBATCH --nodes=1

#SBATCH --partition=batch

#SBATCH --mail-user=eric.tytell@tufts.edu

# Finally, Start the program
module load matlab
matlab -nojvm -nodisplay < run_test_two_muscles.m > test_two_muscles.out