#!/bin/bash 
#SBATCH --job-name=9-admixture
#SBATCH --output=results/slurm_logs/esu/9.admixture.%j.out
#SBATCH --error=results/slurm_logs/esu/9.admixture.%j.err
#SBATCH -t 06:00:00
#SBATCH --partition=shas
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --mail-type=END
#SBATCH  --mail-user=caitlinv@colostate.edu

set -x

module purge
source ~/.bashrc
conda activate admix
cd $1

vcf=results/impute/cawa.sorted.imputedEWCS.miss80.qual30.fst85.filtRda.recode.vcf.gz
outDir=results/esu/admix
name=cawa.imputed.miss80.qual30.fst85.filtRda.run5
num=12345

##plink for input to admixture
##$2 is imputed genotypes
plink --vcf $vcf --recode12 --out $outDir/$name --allow-extra-chr --allow-no-sex --no-sex --no-parents --no-fid --no-pheno --set-missing-var-ids @:#[cawa]\$1,\$2


##Find SNPs in LD
plink --file $outDir/"$name" --out $outDir/LD.snps.exclude.2  --indep-pairwise 200 50 0.5 --allow-extra-chr
##Remove SNPS in LD
plink --file $outDir/"$name" --extract $outDir/LD.snps.exclude.2.prune.in  --recode12 \
	 --out $outDir/"$name".ldremove --allow-extra-chr

##Remove scaffold from map file
sed 's/^scaffold//g' $outDir/"$name".ldremove.map > $outDir/temp.map
mv $outDir/temp.map $outDir/"$name".ldremove.map

##admixture
for K in {1..6}
do
  admixture -j8 -s $num --cv $outDir/"$name".ldremove.ped $K | tee $outDir/$name.${K}.log
done
