#!/bin/bash 
#SBATCH --job-name=imputeAll
#SBATCH --output=results/slurm_logs/impute/impute.%j.out
#SBATCH --error=results/slurm_logs/impute/impute.%j.err
#SBATCH -t 48:00:00
#SBATCH --partition=shas
#SBATCH --qos=long
#SBATCH --ntasks=10
#SBATCH --mem=40G
#SBATCH --mail-type=ALL
#SBATCH  --mail-user=caitlinv@colostate.edu


set -x

module purge
source ~/.bashrc
conda activate beagle

cd $1
outDir="results/impute"
inDir="results/"
outname=cawa.bqsr.miss0.2.qual30.1X.fst85.imputeAll

#Make an exclude list, no single snp scaffolds for beagle
gunzip $inDir/cawa.bqsr.miss0.2.qual30.1X.fst85.recode.vcf.gz
grep -v "^#" $inDir/cawa.bqsr.miss0.2.qual30.1X.fst85.recode.vcf | cut -f 1 | sort | uniq -c | awk '{if ($1 < 10) print $2}'> $outDir/scafExclude.list
grep -w  --file=$outDir/scafExclude.list $inDir/cawa.bqsr.miss0.2.qual30.1X.fst85.recode.vcf| grep -v '^#'  | awk '{print $1":"$2}' > $outDir/excludeMarkers.list


BEAGLE_JAR="/projects/caitlinv@colostate.edu/software/beagle4.1/beagle.27Jan18.7e1.jar"

##impute vcf
java -Xmx40G -jar $BEAGLE_JAR gl=$inDir/cawa.bqsr.miss0.2.qual30.1X.fst85.recode.vcf \
	out=$outDir/$outname.imputed.gl \
	impute=true gp=true \
	nthreads=10 excludemarkers=$outDir/excludeMarkers.list
