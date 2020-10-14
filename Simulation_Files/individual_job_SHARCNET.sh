#!/bin/bash

# # Brendon Phillips
# # PhD Candidate
# # Bauch computational epidemiology research group
# # Department of Applied Mathematics
# # Faculty of Mathematics
# # Universiity of Waterloo

#SBATCH --mem=2G
#SBATCH --output=Jobs/job.%N.%j.out
#SBATCH --ntasks=1

./Skynet_hesitance `cat $1 | head -n $2 | tail -n 1` $3 $4
