#!/bin/bash
##################
#set a job name
#SBATCH --job-name=3b.cov
#################
#a file for job output, you can check job progress
#SBATCH --output=results/slurm_logs/cov/3b.cov.%j.out
#################
# a file for errors from the job
#SBATCH --error=results/slurm_logs/cov/3b.cov.%j.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 48:00:00
#################
#default is one core, should be fine for samtools
#SBATCH --partition=shas
#SBATCH --qos=long
#################
#SBATCH --mem=16G
#################
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
source ~/.bashrc
cd $1

touch results/cov/covSummary.txt

for sample in `ls *.genomecov.txt`
  	do
        awk '{if($3<500) {total+=$3; ++lines}} END {print FILENAME," ",total/lines}' \
        $sample >> results/cov/covSummary.txt
done



