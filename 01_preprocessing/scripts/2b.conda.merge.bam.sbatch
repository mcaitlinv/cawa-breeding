#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=2b.merge.bams
#################
#a file for job output, you can check job progress
#SBATCH --output=results/slurm_logs/bams/2b.merge.bams.%j.out
#################
# a file for errors from the job
#SBATCH --error=results/slurm_logs/bams/2b.merge.bams.%j.err
#################
# a file for errors from the job
#SBATCH --array=1-181
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 01:30:00
#################
#number of nodes
#SBATCH --nodes=1
#SBATCH --ntasks=8
#################
#SBATCH --mem=16G
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
eval $(line_assign.sh $SLURM_ARRAY_TASK_ID samples.txt)
source ~/.bashrc

cd $1
conda activate bwasam
##Get all bams of a single sample name
ls -d results/bams/*"$sample"*.mkdup.bam > results/bams/"$sample".txt

wcBam=$(wc -l < results/bams/"$sample".txt)

if expr "$wcBam" : '[1]'
then
echo "Single bam, skipping merge"
mv  results/bams/"$sample"*.mkdup.bam results/bams/"$sample".mergeMkDup.bam
mv  results/bams/"$sample"*.mkdup.bam.bai results/bams/"$sample".mergeMkDup.bam.bai
else
##merge all bams from file of bams
samtools merge --threads 8 -f results/bams/"$sample".merged.bam -b results/bams/"$sample".txt

##Removing PCR duplicates from BAM files
#sort by name, not position
samtools sort -n -o  results/bams/"$sample".namesort.bam  results/bams/"$sample".merged.bam
rm  results/bams/"$sample".merged.bam

#Add mate score tags for samtools markdup to select best reads
samtools fixmate -m  results/bams/"$sample".namesort.bam  results/bams/"$sample".fixm.bam
rm results/bams/"$sample".namesort.bam

#Sort again by posiition
samtools sort -o  results/bams/"$sample".fixm.sort.bam   results/bams/"$sample".fixm.bam
rm  results/bams/"$sample".fixm.bam

#Markdups
samtools markdup  results/bams/"$sample".fixm.sort.bam  results/bams/"$sample".mergeMkDup.bam
samtools index  results/bams/"$sample".mergeMkDup.bam
rm  results/bams/"$sample".fixm.sort.bam
fi


