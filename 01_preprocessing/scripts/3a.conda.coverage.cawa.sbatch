#!/bin/bash
#SBATCH --job-name=3a.cov
#SBATCH --output=results/slurm_logs/cov/3a.cov.%j.out
#SBATCH --error=results/slurm_logs/cov/3a.cov.%j.err
#SBATCH --time=04:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#SBATCH --mem=16G
#SBATCH --array=1-181
#SBATCH --mail-type=END
#SBATCH --mail-user=caitlinv@colostate.edu

set -x
module purge
eval $(line_assign.sh $SLURM_ARRAY_TASK_ID samples.txt)
source ~/.bashrc
cd $1

##Activate conda
conda activate bedtools

##get depth at each position and save
bedtools genomecov -d -ibam results/bams/"$sample".mergeMkDup.bam  > results/cov/"$sample".genomecov.txt
