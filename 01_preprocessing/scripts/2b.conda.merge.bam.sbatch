#!/bin/bash
#SBATCH --job-name=2b.merge.bams
#SBATCH --output=results/slurm_logs/bams/2b.merge.bams.%j.out
#SBATCH --error=results/slurm_logs/bams/2b.merge.bams.%j.err
#SBATCH --array=1-181
#SBATCH --time=01:30:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#SBATCH --ntasks=8
#SBATCH --mem=16G
#SBATCH --mail-type=END
#SBATCH --mail-user=caitlinv@colostate.edu

set -x
module purge
eval $(line_assign.sh $SLURM_ARRAY_TASK_ID samples.txt)
source ~/.bashrc

cd $1
conda activate bwasam

##Get all bams of a single sample name
ls -d results/bams/*"$sample"*.mkdup.bam > results/bams/"$sample".txt

##Get count of all bams for each sample
wcBam=$(wc -l < results/bams/"$sample".txt)

##if only 1 bam, just rename and move
##else merge multiple bams and move final bam
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
