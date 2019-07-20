#!/bin/bash
#SBATCH -D ./
#SBATCH -J py_samtool
#SBATCH -o get_fragmets4_%A_%a.out
#SBATCH --partition=Lewis
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
##SBATCH --time=2:00:00   #2 hours maybe not enough use 2days
#SBATCH -t 2-00:00
#module add python/python-2.7.14
module add samtools/samtools-1.7 
module add bcftools/bcftools-1.7
module add htslib/htslib-1.7
srun hostname -s | sort -u >slurm.hosts
#echo "ok check this test"
python get_fragments.py 
# > getFrag.log 2>1 &
