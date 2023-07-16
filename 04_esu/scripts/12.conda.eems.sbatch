#!/bin/bash  
#SBATCH --job-name=12.eems
#SBATCH --output=results/slurm_logs/au/12.eems.%j.out
#SBATCH --error=results/slurm_logs/au/12.eems.%j.err
#SBATCH -t 04:00:00
#SBATCH --partition=shas
#SBATCH --nodes=1
#SBATCH --mem=24G
#SBATCH --mail-type=END
#SBATCH  --mail-user=caitlinv@colostate.edu
set -x

module purge
source ~/.bashrc
cd $1

vcf=results/au/cawa.adapt.all.impute.vcf.gz
outDir=results/au/eems
name="cawa.imputed.adapt.all"
num=14944

##plink to get file for bed2diffs
conda activate admix

##plink for input to eems
plink --vcf $vcf --recode --out $outDir/$name --allow-extra-chr --allow-no-sex --no-sex --no-parents --no-fid --no-pheno --set-missing-var-ids @:#[cawa]\$1,\$2

##Dummy code scaffold # as 1, since plink/bed2diffs don't like scaffolds
awk 'BEGIN { FS = "\t"};{ OFS="\t" } ; {print 1, $2, $3, $4, $5}' $outDir/$name.map > $outDir/temp.map
mv $outDir/temp.map $outDir/$name.map
##Get the bed/bim files
plink --file $outDir/$name --out $outDir/$name.$num --make-bed --thin-count $num

##Deactivate conda since eems doesn't play well
conda deactivate admix

##Get genetic dissimilarity matrix
/projects/caitlinv@colostate.edu/software/eems/bed2diffs/src/bed2diffs_v1 --bfile $outDir/$name.$num --nthreads 2

##Params requires names to be the same as in-file
cp $outDir/cawa.ind.coord2 $outDir/$name.$num.coord
cp $outDir/cawa.coord.outer2 $outDir/$name.$num.outer

##make params file
echo "datapath = /scratch/summit/caitlinv@colostate.edu/cawa-breed-wglc/results/au/eems/$name.$num
mcmcpath = /scratch/summit/caitlinv@colostate.edu/cawa-breed-wglc/results/au/eems/data/adapt-all-nInd169-nSites$num-EEMS-nDemes200
nIndiv = 169
nSites = $num
nDemes = 200
diploid = true
numMCMCIter = 10000000
numBurnIter = 10000
numThinIter = 9999" > $outDir/params.$num

##Run eems for snps
/projects/caitlinv@colostate.edu/software/eems/runeems_snps/src/runeems_snps --params $outDir/params.$num --seed 123

