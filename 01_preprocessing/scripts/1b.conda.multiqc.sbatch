#!/bin/bash 
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling  
#################
#set a job name  
#SBATCH --job-name=1b.multiqc
#################  
#a file for job output, you can check job progress
#SBATCH --output=results/slurm_logs/qc/1b.multiqc.%j.out
#################
# a file for errors from the job
#SBATCH --error=results/slurm_logs/qc/1b.multiqc.%j.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 00:30:00
################
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
#SBATCH --partition=shas
#################
#number of nodes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=END
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=caitlinv@colostate.edu
#################
#echo commands to stdout
set -x

module purge
source ~/.bashrc
cd $1

##Activate conda
conda activate multiqc

mkdir -p results/qc/multiqc

multiqc -f -v -o results/qc/multiqc results/qc/fqFASTQC results/qc/trimFASTQC

