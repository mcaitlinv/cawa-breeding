#!/bin/bash
#set a job name
#SBATCH --job-name=5a.mergeVcfs
#SBATCH --output=results/slurm_logs/variants/5a.mergeVcfs.%j.out
#SBATCH --error=results/slurm_logs/variants/5a.mergeVcfs.%j.err
################
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#################
# Note: 4.84G/core or task
#################
#SBATCH --mem=32G
#################
set -x

module purge
source ~/.bashrc
cd $1

conda activate bcftools

##Make temp dir, bcftools was not finding default tmp
mkdir -p tmp

# list file with the 535 2mb intervals
vcfs=intervals2mb.vcfs.list
bcfs=intervals2mb.bcfs.list

##out names for merged files
outname=cawa.unfiltered

##index bcfs
while read a; do bcftools index $a; done < intervals2mb.bcfs.list
##bcftools concat for bcfs
bcftools concat -a  --threads 8 -f ${bcfs} -Oz -o results/variants/"$outname".bcf.gz
bcftools sort -m 30G -T tmp -Oz -o results/variants/"$outname".sort.bcf.gz  results/variants/"$outname".bcf.gz

##No need to index, gatk already made it
##bcftools concat for bcfs
bcftools concat -a  --threads 8 -f ${vcfs} -Oz -o results/variants/"$outname".vcf.gz
bcftools sort -m 30G -T tmp -Oz -o results/variants/"$outname".sort.vcf.gz  results/variants/"$outname".vcf.gz



