#!/bin/bash
#SBATCH --job-name=7a.mergeVcfs
#SBATCH --output=results/slurm_logs/variants/7a.mergeVcfs.%j.out
#SBATCH --error=results/slurm_logs/variants/7a.mergeVcfs.%j.err
#SBATCH --time=12:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#SBATCH --mem=32G

set -x
module purge
source ~/.bashrc
cd $1
conda activate bcftools

##Make temp dir, bcftools was not finding default tmp
mkdir -p tmp

# list file with the 535 2mb intervals
vcfs=intervals2mb.bqsr.vcfs.list

#out dir
outDir=results/variants2
##out names for merged files
outname=cawa.bqsr.unfiltered

##No need to index, gatk already made it
##bcftools concat for vcfs
bcftools concat -a  --threads 8 -f $vcfs -Oz -o $outDir/"$outname".vcf.gz
bcftools sort -m 30G -T tmp -Oz -o $outDir/"$outname".sort.vcf.gz  $outDir/"$outname".vcf.gz
