#!/bin/bash
#set a job name
#SBATCH --job-name=get.mergeVcfs
#SBATCH --output=results/slurm_logs/postprocess/get.mergeVcfs.%j.out
#SBATCH --error=results/slurm_logs/postprocess/get.mergeVcfs.%j.err
################
#SBATCH --time=04:00:00
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
vcfs=fst2mb.vcfs.list

##out names for merged files
outname=fst.unfiltered

##No need to index, gatk already made it
##bcftools concat for bcfs
bcftools concat -a  --threads 8 -f ${vcfs} -Oz -o results/fstbams/vcfs/"$outname".vcf.gz
bcftools sort -m 30G -T tmp -Oz -o results/fstbams/vcfs/"$outname".sort.vcf.gz  results/fstbams/vcfs/"$outname".vcf.gz



