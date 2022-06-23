#!/bin/bash 
#SBATCH --job-name=fstRemove
#SBATCH --output=results/slurm_logs/postprocess/fstRemove.%j.out
#SBATCH --error=results/slurm_logs/postprocess/fstRemove.%j.err
#SBATCH --partition=shas
#SBATCH -t 03:00:00
#SBATCH --nodes=1
#SBATCH --mem=24G
#SBATCH --mail-type=END
#SBATCH  --mail-user=caitlinv@colostate.edu

set -x

module purge
source ~/.bashrc
cd $1
conda activate vcftools

##directory
outDir="results/fstbams/vcfs"
vcf=fst.unfiltered.sort.vcf.gz

##lists of sample names in each "population"
##really just same samples on different platforms
touch hiseq.samples.list;for f in results/fstbams/*plate{1..2}*RG.bam; do echo `basename $f` | cut -f1,2 -d '.' >> hiseq.samples.list; done
touch novaseq.samples.list; for f in results/fstbams/*plate3*RG.bam; do echo `basename $f` | cut -f1,2 -d '.' >> novaseq.samples.list; done

pop1=novaseq.samples.list
pop2=hiseq.samples.list

##get fst between platforms
vcftools --gzvcf $outDir/$vcf --weir-fst-pop $pop1 --weir-fst-pop $pop2 --out $outDir/fst.outliers.hiseq.novaseq

##Change -nan values to 0
sed -i 's/-nan/0/g' $outDir/fst.outliers.hiseq.novaseq.weir.fst

# ##Get 90-99th percentiles of fst values
awk '{if (NR>1) {print $3}}' $outDir/fst.outliers.hiseq.novaseq.weir.fst > $outDir/fst.txt
pc80=$(cat $outDir/fst.txt|(percentile=80; (sort -n;echo)|nl -ba -v0|tac|(read count;cut=$(((count * percentile + 99) / 100)); tac|sed -n "${cut}s/.*\t//p")))
pc85=$(cat $outDir/fst.txt|(percentile=85; (sort -n;echo)|nl -ba -v0|tac|(read count;cut=$(((count * percentile + 99) / 100)); tac|sed -n "${cut}s/.*\t//p")))
pc90=$(cat $outDir/fst.txt|(percentile=90; (sort -n;echo)|nl -ba -v0|tac|(read count;cut=$(((count * percentile + 99) / 100)); tac|sed -n "${cut}s/.*\t//p")))
pc95=$(cat $outDir/fst.txt|(percentile=95; (sort -n;echo)|nl -ba -v0|tac|(read count;cut=$(((count * percentile + 99) / 100)); tac|sed -n "${cut}s/.*\t//p")))
pc99=$(cat $outDir/fst.txt|(percentile=99; (sort -n;echo)|nl -ba -v0|tac|(read count;cut=$(((count * percentile + 99) / 100)); tac|sed -n "${cut}s/.*\t//p")))

# ##pull out variants above 90-99th percentiles
awk -v var="$pc80" 'NR>1 && $3 > var' $outDir/fst.outliers.hiseq.novaseq.weir.fst > $outDir/fstfilter.80.txt
awk -v var="$pc85" 'NR>1 && $3 > var' $outDir/fst.outliers.hiseq.novaseq.weir.fst > $outDir/fstfilter.85.txt
awk -v var="$pc90" 'NR>1 && $3 > var' $outDir/fst.outliers.hiseq.novaseq.weir.fst > $outDir/fstfilter.90.txt
awk -v var="$pc95" 'NR>1 && $3 > var' $outDir/fst.outliers.hiseq.novaseq.weir.fst > $outDir/fstfilter.95.txt
awk -v var="$pc99" 'NR>1 && $3 > var' $outDir/fst.outliers.hiseq.novaseq.weir.fst > $outDir/fstfilter.99.txt
