#!/bin/bash 
#SBATCH --job-name=8c.impute
#SBATCH --output=results/slurm_logs/impute/8c.impute.%j.out
#SBATCH --error=results/slurm_logs/impute/8c.impute.%j.err
#SBATCH -t 2:00:00
#SBATCH --partition=shas
#SBATCH --ntasks=8
#SBATCH --mail-type=END
#SBATCH  --mail-user=caitlinv@colostate.edu

set -x

module purge
source ~/.bashrc
conda activate bcftools

cd $1
outDir="results/impute"


##Merge vcfs back together
touch $outDir/bgzip.list
ls $outDir/*.gl.vcf.gz > $outDir/bgzip.list
while read a; do
bcftools index $a; done < $outDir/bgzip.list

##this merges in sample order given by vcfs, not sorted ordered
bcftools merge -l $outDir/bgzip.list --threads 8 --force-samples -o $outDir/cawa.imputedEWCS.miss80.qual30.vcf.gz

##So sort the samples
bcftools query -l $outDir/cawa.imputedEWCS.miss80.qual30.vcf.gz | sort > $outDir/sort.samples.txt
bcftools view -S $outDir/sort.samples.txt $outDir/cawa.imputedEWCS.miss80.qual30.vcf.gz > $outDir/cawa.sorted.imputedEWCS.miss80.qual30.vcf.gz


