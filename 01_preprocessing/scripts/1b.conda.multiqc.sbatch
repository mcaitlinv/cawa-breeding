#!/bin/bash 
#SBATCH --job-name=1b.multiqc
#SBATCH --output=results/slurm_logs/qc/1b.multiqc.%j.out
#SBATCH --error=results/slurm_logs/qc/1b.multiqc.%j.err
#SBATCH --time=00:30:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mail-type=END
#SBATCH --mail-user=caitlinv@colostate.edu

set -x
module purge
source ~/.bashrc
cd $1

##Activate conda
conda activate multiqc

mkdir -p results/qc/multiqc

##Get multiqc file
multiqc -f -v -o results/qc/multiqc results/qc/fqFASTQC results/qc/trimFASTQC
