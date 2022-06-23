#!/bin/bash
#SBATCH --job-name=2a.bwa
#SBATCH --output=results/slurm_logs/map/2a.bwa.%j.out
#SBATCH --error=results/slurm_logs/map/2a.bwa.%j.err
#SBATCH --time=03:00:00
#SBATCH --array=1-299
#SBATCH --qos=normal
#SBATCH --partition=shas
#SBATCH --ntasks-per-node 8
#SBATCH --mem=16G
#SBATCH --mail-type=END
#SBATCH --mail-user=caitlinv@colostate.edu

set -x
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
