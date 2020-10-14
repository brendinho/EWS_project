#!/bin/bash

# # Brendon Phillips
# # PhD Candidate
# # Bauch computational epidemiology research group
# # Department of Applied Mathematics
# # Faculty of Mathematics
# # Universiity of Waterloo

#SBATCH --job-name="dopeshit"
#SBATCH --partition=hagrid_long
#SBATCH --account=hagrid
#SBATCH --mem=2G
#SBATCH --mail-user=b2philli@uwaterloo.ca
# #SBATCH -- mail-type=FAIL
# #SBATCH --workdir="/work/b2philli"
#SBATCH --output=Jobs/host_name-%j.out
#SBATCH --error=Jobs/host_name-%j.err
#SBATCH --ntasks=1

./Skynet_hesitance `cat $1 | head -n $2 | tail -n 1` $3 $4
