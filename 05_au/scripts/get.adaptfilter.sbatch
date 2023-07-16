#!/bin/bash
#set a job name
#SBATCH --job-name=adaptfiltervcf
#SBATCH --output=results/slurm_logs/au/adaptfiltervcf.%j.out
#SBATCH --error=results/slurm_logs/au/adaptfiltervcf.%j.err
#SBATCH --time=00:30:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#SBATCH --nodes=1
#SBATCH --mem=32G
set -x

module purge
source ~/.bashrc
cd $1


##Make temp dir, bcftools was not finding default tmp
mkdir -p tmp

##the vcf to filter
impute=results/impute/cawa.sorted.imputedEWCS.miss80.qual30.fst85.filtRda.recode.vcf.gz

##out directory
outDir="results/au"

##filter vcfs
outname=cawa.adapt.all

##adaptive vars file
file=resources/adaptive.lfmm.01.rda.SD3.snps.txt

conda activate bcftools

bcftools view --threads 8 --targets-file $file -Oz -o $outDir/$outname.vcf.gz $impute

rm -r tmp
