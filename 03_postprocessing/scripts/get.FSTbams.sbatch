#!/bin/bash 
#SBATCH --job-name=getFSTBams
#SBATCH --output=results/slurm_logs/postprocess/getFSTBams.%j.out
#SBATCH --error=results/slurm_logs/postprocess/getFSTBams.%j.err
################
#SBATCH --partition=shas
#SBATCH -t 04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node 4
#SBATCH --mem=16G
#SBATCH --mail-type=END
#SBATCH  --mail-user=caitlinv@colostate.edu

set -x

module purge
source ~/.bashrc
cd $1


bamDir="results/bqsrBams"
fstDir="results/fstbams"

mkdir -p results/fstbams

while read a b; do 
   cp $bamDir/"$b" $fstDir/"$b"
 done  < fstRemoveBams.txt

conda activate bwasam
##find read groups
cd $fstDir

for f in *.recal.bam;do \
 samtools split $f -f %!
 done

rm *.recal.bam

for f in *; do \
  mv $f $f.recal.bam ;
  samtools index $f.recal.bam
done
 
cd $1 
conda activate picard

for f in $fstDir/*recal.bam; do 
  sample=$(echo `basename $f` | cut -f 1,2 -d'.') ;
  id=$(echo `basename $f`| cut -f 1,2,3,4 -d'.') ;
  plate=$(echo `basename $f` | cut -f 2 -d'.') ;
  lane=$(echo `basename $f` | cut -f 3,4 -d'.');
  picard AddOrReplaceReadGroups -I $f \
  -O $fstDir/`basename $f .recal.bam`.newRG.bam \
  -RGID $id -RGLB illumina.$lane -RGPL $plate \
  -RGSM $sample -RGPU $sample --CREATE_INDEX true
done


