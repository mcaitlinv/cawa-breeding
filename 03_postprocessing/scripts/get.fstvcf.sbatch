#!/bin/bash
#set a job name
#SBATCH --job-name=get.fstfiltervcf
#SBATCH --output=results/slurm_logs/postprocess/get.fstfiltervcf.%j.out
#SBATCH --error=results/slurm_logs/postprocess/get.fstfiltervcf.%j.err
################
#SBATCH --time=03:00:00
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

##the vcf to filter
vcf=results/variants/cawa.bqsr.unfiltered.sort.vcf.gz

##out directory
outDir="results/fstbams/vcfs"

##filter vcfs
outname=cawa.bqsr.miss0.2.qual30.1X
bcftools view --threads 8 -S samples1X.list -m 2 -M 2 --min-af 0.05 --max-af 0.95 -i 'F_MISSING < 0.2 && QUAL>30' -Oz -o $outDir/$outname.vcf.gz $vcf

conda activate vcftools
vcftools --gzvcf $outDir/$outname.vcf.gz --recode --keep-INFO-all --exclude-positions $outDir/fstfilter.80.txt --out $outDir/$outname.fst80
vcftools --gzvcf $outDir/$outname.vcf.gz --recode --keep-INFO-all --exclude-positions $outDir/fstfilter.85.txt --out $outDir/$outname.fst85
vcftools --gzvcf $outDir/$outname.vcf.gz --recode --keep-INFO-all --exclude-positions $outDir/fstfilter.90.txt --out $outDir/$outname.fst90
vcftools --gzvcf $outDir/$outname.vcf.gz --recode --keep-INFO-all --exclude-positions $outDir/fstfilter.95.txt --out $outDir/$outname.fst95
vcftools --gzvcf $outDir/$outname.vcf.gz --recode --keep-INFO-all --exclude-positions $outDir/fstfilter.99.txt --out $outDir/$outname.fst99

conda activate bcftools
##get allele depths for SRS PCA analysis
bcftools query -f '%CHROM\t%POS[\t%AD]\n' $outDir/$outname.vcf.gz > $outDir/$outname.AD.txt
bcftools query -f '%CHROM\t%POS[\t%AD]\n' $outDir/$outname.fst80.recode.vcf > $outDir/$outname.fst80.AD.txt
bcftools query -f '%CHROM\t%POS[\t%AD]\n' $outDir/$outname.fst85.recode.vcf > $outDir/$outname.fst85.AD.txt
bcftools query -f '%CHROM\t%POS[\t%AD]\n' $outDir/$outname.fst90.recode.vcf > $outDir/$outname.fst90.AD.txt
bcftools query -f '%CHROM\t%POS[\t%AD]\n' $outDir/$outname.fst95.recode.vcf > $outDir/$outname.fst95.AD.txt
bcftools query -f '%CHROM\t%POS[\t%AD]\n' $outDir/$outname.fst99.recode.vcf > $outDir/$outname.fst99.AD.txt

rm -r tmp
