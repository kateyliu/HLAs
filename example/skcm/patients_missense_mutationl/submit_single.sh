#!/bin/sh
#SBATCH --partition Lewis
#SBATCH --ntasks 1
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 1G
#SBATCH -t 2-00:00
#SBATCH -o step5_%A_%a.out


source ~/conda.bashrc
source activate py2
#add the py env

python get_amp_single.py $1
