#!/bin/bash
#set a job name
#SBATCH --job-name=4a.mpileup
#SBATCH --output=results/slurm_logs/variants/4a.mpileup.%A.%a.out
#SBATCH --error=results/slurm_logs/variants/4a.mpileup.%A.%a.err
################
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#################
# Note: 4.84G/core or task
#################
#SBATCH --mem=16G
#SBATCH --array=1-535
#################
set -x

module purge
source ~/.bashrc
conda activate bcftools
cd $1

##Get interval file name
intervals=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' intervals2mb.txt)
##Get interval set number
intervalID=$(echo ${intervals} | cut -f2 -d_ | cut -f1 -d.)

##set out name
mkdir -p results/variants/mpileup
outname=$(printf 'results/variants/mpileup/cawa2mb.interval.%s.samtools.bcf'  $intervalID)

##Get list of bams
ls results/bams/*mergeMkDup.bam > cawa.all.bam.list
bamlist=cawa.all.bam.list
reference="/projects/caitlinv@colostate.edu/genomes/cardellina_canadensis_pseudohap_v1.fasta"

bcftools mpileup -Ou -f $reference -R ${intervals} --annotate FORMAT/AD,FORMAT/DP -b ${bamlist} | bcftools call -mv -Ob -o ${outname}

