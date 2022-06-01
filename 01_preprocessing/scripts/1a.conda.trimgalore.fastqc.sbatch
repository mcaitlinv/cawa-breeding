#!/bin/bash 
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling  
#################
#set a job name  
#SBATCH --job-name=1a.trim.fastqc
#################  
#a file for job output, you can check job progress
#SBATCH --output=results/slurm_logs/trim/1a.trim.fastqc.%j.out
#################
# a file for errors from the job
#SBATCH --error=results/slurm_logs/trim/1a.trim.fastqc.%j.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 01:00:00
#################
#SBATCH --array=1-299
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
#SBATCH --partition=shas
#################
#number of nodes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
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
eval $(line_assign.sh $SLURM_ARRAY_TASK_ID numbered-units.txt) 
cd $1

##Activate conda
conda activate trimQC

mkdir -p results/slurm_logs/trim/
mkdir -p results/trim/
mkdir -p results/qc/fqFASTQC
mkdir -p results/logs/qc/fqFASTQC
mkdir -p results/qc/trimFASTQC
mkdir -p results/logs/qc/trimFASTQC

T1=`basename $fq1 .fq.gz`_val_1.fq.gz
T2=`basename $fq2 .fq.gz`_val_2.fq.gz
O1=`basename $T1 _val_1.fq.gz`.trim.fq.gz
O2=`basename $T2 _val_2.fq.gz`.trim.fq.gz

fastqc -o results/qc/fqFASTQC --noextract fastqs/$fq1 fastqs/$fq2 > results/logs/qc/fqFASTQC/fastqc."$fq1".log 2>&1



###Trim low quality fragments and adapters
trim_galore -q 5 --paired --cores 8 --output_dir results/trim fastqs/$fq1 fastqs/$fq2


##Fastp for trimming the poly-g tails
/projects/caitlinv@colostate.edu/software/fastp --in1  results/trim/$T1 --out1 results/trim/$O1 \
	--in2 results/trim/$T2 --out2 results/trim/$O2  --thread 10 \
	-L -A --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20

fastqc -o results/qc/trimFASTQC --noextract results/trim/$O1 results/trim/$O2 > results/logs/qc/trimFASTQC/fastqc."$O1".log 2>&1
