#!/bin/bash
#SBATCH --job-name=6a.bqsr
#SBATCH --output=results/slurm_logs/bqsr/6a.bqsr.%j.out
#SBATCH --error=results/slurm_logs/bqsr/6a.bqsr.%j.err
#SBATCH --time=02:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#SBATCH --array=1-181
#SBATCH --mem=32G

set -x
module purge
source ~/.bashrc
eval $(line_assign.sh $SLURM_ARRAY_TASK_ID samples.txt)
cd $1
conda activate gatk4.2.5.0

##Make dir to put tables/bams
mkdir -p results/bqsr

##these are all the input/ouput names and directories
inDir=results/bams
bam=$sample.mergeMkDup.bam
outtable=$sample.table.out
outDir=results/bqsr

##ref genome
reference="/projects/caitlinv@colostate.edu/genomes/cardellina_canadensis_pseudohap_v1.fasta"

##get recal tables for each sample
##known sites is intersection of samtools/gatk calls and hard filtered
gatk --java-options "-Xmx30g -XX:+UseParallelGC -XX:ParallelGCThreads=8" BaseRecalibrator \
-I $inDir/$bam --known-sites results/variants/isec/0002.vcf.gz \
-O $outDir/$outtable -R $reference

##recal each sample bam using the recal table
gatk --java-options "-Xmx30g -XX:+UseParallelGC -XX:ParallelGCThreads=8" ApplyBQSR \
-R $reference -I $inDir/$bam \
--bqsr-recal-file $outDir/$outtable -O $outDir/$sample.recal.bam

##index new bams
conda activate bwasam
samtools index $outDir/$sample.recal.bam

##Get what recal table looks like for new bams
gatk --java-options "-Xmx30g -XX:+UseParallelGC -XX:ParallelGCThreads=8" BaseRecalibrator \
-I $outDir/$sample.recal.bam --known-sites results/variants/isec/0002.vcf.gz \
-O $outDir/$outtable.after -R $reference

##Get the graphs of what bqsr did to qualities
gatk --java-options "-Xmx30g -XX:+UseParallelGC -XX:ParallelGCThreads=8" AnalyzeCovariates \
     -R $reference \
     -before $outDir/$outtable \
     -after $outDir/$outtable.after \
     -plots $outDir/$sample.analyzeCovariates.pdf
