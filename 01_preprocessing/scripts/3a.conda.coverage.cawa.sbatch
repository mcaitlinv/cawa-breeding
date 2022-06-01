#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=3a.cov
#################
#a file for job output, you can check job progress
#SBATCH --output=results/slurm_logs/cov/3a.cov.%j.out
#################
# a file for errors from the job
#SBATCH --error=results/slurm_logs/cov/3a.cov.%j.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 04:00:00
#################
#default is one core, should be fine for samtools
#################
#SBATCH --mem=16G
#################
#SBATCH --array=1-181
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=END
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=caitlinv@colostate.edu
##################
#echo commands to stdout
set -x

module purge
eval $(line_assign.sh $SLURM_ARRAY_TASK_ID samples.txt)
source ~/.bashrc
cd $1

##Activate conda
conda activate bedtools

bedtools genomecov -d -ibam results/bams/"$sample".mergeMkDup.bam  > results/cov/"$sample".genomecov.txt
