#!/bin/bash
#SBATCH --job-name=7b.filterVcfs
#SBATCH --output=results/slurm_logs/variants/7b.filterVcfs.%j.out
#SBATCH --error=results/slurm_logs/variants/7b.filterVcfs.%j.err
#SBATCH --time=4:00:00
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

##the vcf to filter
vcf=cawa.bqsr.unfiltered.sort.vcf.gz

##out directory
outDir="results/variants"

##index the vcfs
bcftools index $outDir/$vcf

##Get unfiltered stats
bcftools stats $outDir/$vcf > $outDir/$vcf.stats

##filter vcfs
outname=cawa.bqsr.miss0.2.qual30
bcftools view --threads 8 -m 2 -M 2 -S samples1X.list --min-af 0.05 --max-af 0.95 -i 'F_MISSING < 0.2 && QUAL>30' -Oz -o $outDir/$outname.1X.vcf.gz $outDir/$vcf

##index filtered
bcftools index $outDir/"$outname".1X.vcf.gz

##stats on filtered
bcftools stats $outDir/"$outname".vcf.gz > $outDir/"$outname".1X.vcf.stats

##get allele depths for SRS PCA analysis
bcftools query -f '%CHROM\t%POS[\t%AD]\n' $outDir/"$outname".1X.vcf.gz > $outDir/"$outname".1X.AD.txt
