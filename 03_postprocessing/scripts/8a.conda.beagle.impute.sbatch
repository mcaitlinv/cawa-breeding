#!/bin/bash 
#SBATCH --job-name=8a.impute
#SBATCH --output=results/slurm_logs/impute/8a.impute.%j.out
#SBATCH --error=results/slurm_logs/impute/8a.impute.%j.err
#SBATCH -t 12:00:00
#SBATCH --partition=shas
#SBATCH --ntasks=2
#SBATCH --mem=8G
#SBATCH --mail-type=END
#SBATCH  --mail-user=caitlinv@colostate.edu

set -x

module purge
source ~/.bashrc
conda activate bcftools

cd $1

##outdirectory
mkdir -p results/impute
inDir="results/variants"
outDir="results/impute"

##Split entire vcf into population groups
bcftools view -S east.txt --output $outDir/east.filtered.vcf.gz $inDir/cawa.bqsr.miss0.2.qual30.1X.fst85.recode.vcf
bcftools view -S north.txt --output $outDir/north.filtered.vcf.gz $inDir/cawa.bqsr.miss0.2.qual30.1X.fst85.recode.vcf
bcftools view -S south.txt --output $outDir/south.filtered.vcf.gz $inDir/cawa.bqsr.miss0.2.qual30.1X.fst85.recode.vcf

#Make an exclude list, no single snp scaffolds for beagle
grep -v "^#" $inDir/cawa.bqsr.miss0.2.qual30.1X.fst85.recode.vcf | cut -f 1 | sort | uniq -c | awk '{if ($1 < 10) print $2}'> $outDir/scafExclude.list
grep -w  --file=$outDir/scafExclude.list $inDir/cawa.bqsr.miss0.2.qual30.1X.fst85.recode.vcf| grep -v '^#'  | awk '{print $1":"$2}' > $outDir/excludeMarkers.list
#bcftools view --output $inDir/cawa.bqsr.miss0.2.qual30.1X.fst85.vcf.gz $inDir/cawa.bqsr.miss0.2.qual30.1X.fst85.recode.vcf

##Submit jobs to impute missing data
jid1=$(sbatch /projects/caitlinv@colostate.edu/lcwg-scripts/8b.conda.beagle.impute.sbatch $PWD $outDir/east.filtered.vcf.gz | cut -f4 -d" ")
jid2=$(sbatch /projects/caitlinv@colostate.edu/lcwg-scripts/8b.conda.beagle.impute.sbatch $PWD $outDir/north.filtered.vcf.gz | cut -f4 -d" ")
jid3=$(sbatch /projects/caitlinv@colostate.edu/lcwg-scripts/8b.conda.beagle.impute.sbatch $PWD $outDir/south.filtered.vcf.gz | cut -f4 -d" ")

# jid4=$(sbatch --dependency=afterok:$jid1:$jid2:$jid3 \
# 	/projects/caitlinv@colostate.edu/lcwg-scripts/8c.conda.beagle.impute.sbatch $PWD)

