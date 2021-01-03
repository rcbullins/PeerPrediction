#!/bin/bash
#SBATCH --account=peerprediction
#SBATCH --job-name=asmb4
#SBATCH --output=/home/<VT_PID>/Results/Out/asmb4.out
#SBATCH --error=/home/<VT_PID>/Results/Error/asmb4.err
#SBATCH --mail-user <VT_PID>@vt.edu
#SBATCH --mail-type=end
#SBATCH --mem=100g
#SBATCH --partition=normal_q
#SBATCH --time=144:00:00
module load matlab/R2018a
matlab -nodisplay -nosplash -nodesktop -r "run('/home/<VT_PID>/Code/Assembly/asmb_part4.m');exit;"
