#!/bin/bash
#SBATCH --job-name=5b.isecVcfs
#SBATCH --output=results/slurm_logs/variants/5b.isecVcfs.%j.out
#SBATCH --error=results/slurm_logs/variants/5b.isecVcfs.%j.err
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

##the 2 vcfs to filter then intersect
vcf=temp.sort.vcf.gz
bcf=temp.sort.bcf.gz

##out directory
outDir="results/variants"

##index the vcfs 
bcftools index $outDir/$vcf
bcftools index $outDir/$bcf

##Get unfiltered stats
bcftools stats $outDir/$vcf > $outDir/$vcf.stats
bcftools stats $outDir/$bcf > $outDir/$bcf.stats

##filter vcfs
outname=cawa.snps.miss0.1.qual30
bcftools view --threads 8 -m 2 -M 2 --min-af 0.05 --max-af 0.95 -i 'F_MISSING < 0.1 && QUAL>30' -Oz -o $outDir/$outname.vcf.gz $outDir/$vcf
bcftools view --threads 8 -m 2 -M 2 --min-af 0.05 --max-af 0.95 -i 'F_MISSING < 0.1 && QUAL>30' -Oz -o $outDir/$outname.bcf.gz $outDir/$bcf

##index filtered
bcftools index $outDir/"$outname".vcf.gz
bcftools index $outDir/"$outname".bcf.gz

##stats on filtered
bcftools stats $outDir/"$outname".vcf.gz > $outDir/"$outname".vcf.stats
bcftools stats $outDir/"$outname".bcf.gz > $outDir/"$outname".bcf.stats

##intersect variant sets, get stats on intersection
mkdir -p $outDir/isec
bcftools isec $outDir/"$outname".vcf.gz $outDir/"$outname".bcf.gz -p $outDir/isec -Oz
bcftools index $outDir/isec/*
bcftools stats $outDir/isec/0002.vcf.gz > $outDir/isec/0002.vcf.stats
