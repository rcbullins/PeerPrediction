#!/bin/bash
#SBATCH --account=peerprediction
#SBATCH --job-name= asmb2
#SBATCH --output=/home/<VT_PID>/Results/Out/asmb2.out
#SBATCH --error=/home/<VT_PID>/Results/Error/asmb2.err
#SBATCH --mail-user <VT_PID>@vt.edu
#SBATCH --mail-type=end
#SBATCH --mem=100g
#SBATCH --partition=normal_q
#SBATCH --time=144:00:00
module load matlab/R2018a
matlab -nodisplay -nosplash -nodesktop -r "run('/home/<VT_PID>/Code/Assembly/asmb_part2.m');exit;"
