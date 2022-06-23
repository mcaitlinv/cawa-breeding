#!/bin/bash
#SBATCH --job-name=3b.cov
#SBATCH --output=results/slurm_logs/cov/3b.cov.%j.out
#SBATCH --error=results/slurm_logs/cov/3b.cov.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=shas
#SBATCH --qos=long
#SBATCH --mem=16G
#SBATCH --mail-type=END
#SBATCH --mail-user=caitlinv@colostate.edu

set -x
module purge
source ~/.bashrc
cd $1

touch results/cov/covSummary.txt

##get ave depth from position depths and then remove position file
for sample in `ls results/cov/*.genomecov.txt`
  	do
        awk '{if($3<500) {total+=$3; ++lines}} END {print FILENAME," ",total/lines}' \
        $sample >> results/cov/covSummary.txt
done
