#!/bin/bash 
#SBATCH --job-name=adapt.highgrade
#SBATCH --output=results/slurm_logs/au/adapt.highgrade.%j.out
#SBATCH --error=results/slurm_logs/au/adapt.highgrad.%j.err
#SBATCH -t 02:00:00
#SBATCH --partition=shas
#SBATCH --ntasks=4
#SBATCH --mail-type=END
#SBATCH  --mail-user=caitlinv@colostate.edu
set -x

module purge
source ~/.bashrc
conda activate admix
cd $1

vcf=results/au/cawa.adapt.all.vcf.gz
outDir=results/au/admix/highgrade
num=1

for f in {1..5}
do
name=cawa.adapt.highgrade."$f".100
##plink for input to admixture
##$2 is imputed genotypes
plink --vcf $vcf --recode12 --out $outDir/$name --thin-indiv-count 100 --allow-extra-chr --allow-no-sex --no-sex --no-parents --no-fid --no-pheno --set-missing-var-ids @:#[cawa]\$1,\$2

##Remove scaffold from map file
sed 's/^scaffold//g' $outDir/"$name".map > $outDir/temp.map
mv $outDir/temp.map $outDir/"$name".map

##admixture
for K in {1..8}
do
  admixture -j4 -s $num --cv $outDir/"$name".ped $K | tee $outDir/$name.${K}.log

done

done
mv "$name"*.P $outDir/
mv "$name"*.Q $outDir/

grep "CV" $outDir/impute.$name.*.log | awk '{print $3,$4}' | cut -c 4,7-20 > $outDir/$name.cv.error

