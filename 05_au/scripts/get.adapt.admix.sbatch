#!/bin/bash 
#SBATCH --job-name=adapt.admix
#SBATCH --output=results/slurm_logs/au/adapt.admix.%j.out
#SBATCH --error=results/slurm_logs/au/adapt.admix.%j.err
#SBATCH -t 00:15:00
#SBATCH --partition=shas
#SBATCH --ntasks=4
#SBATCH --mail-type=END
#SBATCH  --mail-user=caitlinv@colostate.edu

set -x

module purge
source ~/.bashrc
conda activate admix
cd $1

vcf1=results/au/cawa.adapt.intersect.impute.vcf.gz
vcf2=results/au/cawa.adapt.all.impute.vcf.gz
outDir=results/au/admix
name1=cawa.adapt.intersect.run1
name2=cawa.adapt.all.run1
num=1

##plink for input to admixture
##$2 is imputed genotypes
plink --vcf $vcf1 --recode12 --out $outDir/$name1 --allow-extra-chr --allow-no-sex --no-sex --no-parents --no-fid --no-pheno --set-missing-var-ids @:#[cawa]\$1,\$2
plink --vcf $vcf2 --recode12 --out $outDir/$name2 --allow-extra-chr --allow-no-sex --no-sex --no-parents --no-fid --no-pheno --set-missing-var-ids @:#[cawa]\$1,\$2


##Remove scaffold from map file
sed 's/^scaffold//g' $outDir/"$name1".map > $outDir/temp1.map
mv $outDir/temp1.map $outDir/"$name1".map
sed 's/^scaffold//g' $outDir/"$name2".map > $outDir/temp2.map
mv $outDir/temp2.map $outDir/"$name2".map

##admixture
for K in {1..8}
do
  admixture -j4 -s $num --cv $outDir/"$name1".ped $K | tee $outDir/impute.$name1.${K}.log
  admixture -j4 -s $num --cv $outDir/"$name2".ped $K | tee $outDir/impute.$name2.${K}.log

done

mv "$name1"*.P $outDir/
mv "$name2"*.P $outDir/

mv "$name1"*.Q $outDir/
mv "$name2"*.Q $outDir/

grep "CV" $outDir/impute.$name1.*.log | awk '{print $3,$4}' | cut -c 4,7-20 > $outDir/$name1.cv.error
grep "CV" $outDir/impute.$name2.*.log | awk '{print $3,$4}' | cut -c 4,7-20 > $outDir/$name2.cv.error


