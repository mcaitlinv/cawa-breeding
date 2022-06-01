#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=2a.bwa
#################
#a file for job output, you can check job progress for each sample
#SBATCH --output=results/slurm_logs/map/2a.bwa.%j.out
#################
# a file for errors from the job for each sample
#SBATCH --error=results/slurm_logs/map/2a.bwa.%j.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 03:00:00
#################
#SBATCH --array=1-299
#################
#number of nodes
#SBATCH -N 1
#SBATCH --ntasks-per-node 8
#################
#SBATCH --mem=16G
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=END
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=caitlinv@colostate.edu
#################
#now run normal batch commands
##################
#echo commands to stdout
set -x
##################

module purge 
source ~/.bashrc
eval $(line_assign.sh $SLURM_ARRAY_TASK_ID numbered-units.txt)

cd $1
mkdir -p results/bams

##Variables: trimmed fq names, genome directory, new bam names
R1=`basename $fq1 .fq.gz`.trim.fq.gz
R2=`basename $fq2 .fq.gz`.trim.fq.gz
genomeDIR="/projects/caitlinv@colostate.edu/genomes/cardellina_canadensis_pseudohap_v1.fasta"
ID="$sample"."$library"."$flowcell"."$lane"

##Align each sample to genome. Note that genome reference must already be built through bwa
conda activate bwasam

##map paired reads, trimmed reads
bwa mem  -t 8 $genomeDIR results/trim/$R1 results/trim/$R2 > results/bams/"$ID".sam

#########sort, add read group information and index it#########
samtools sort -o results/bams/"$ID".bam results/bams/"$ID".sam

##Add read groups, and adjust them accordingling
picard AddOrReplaceReadGroups -INPUT results/bams/"$ID".bam -RGID "$ID" -RGLB "$library" \
	-RGPL "$platform"."$flowcell"."$lane" -RGPU "$library"."$sample" -RGSM "$sample" \
	-OUTPUT results/bams/"$ID".RG.bam -VALIDATION_STRINGENCY SILENT
samtools index results/bams/"$ID".RG.bam

##Removing PCR duplicates from BAM files
#sort by name, not position
samtools sort -n -o results/bams/"$ID".namesort.bam results/bams/"$ID".RG.bam

#Add mate score tags for samtools markdup to select best reads
samtools fixmate -m results/bams/"$ID".namesort.bam results/bams/"$ID".fixm.bam

#Sort again by posiition
samtools sort -o results/bams/"$ID".fixm.sort.bam  results/bams/"$ID".fixm.bam

#Markdups
samtools markdup results/bams/"$ID".fixm.sort.bam results/bams/"$ID".mkdup.bam
samtools index results/bams/"$ID".mkdup.bam
