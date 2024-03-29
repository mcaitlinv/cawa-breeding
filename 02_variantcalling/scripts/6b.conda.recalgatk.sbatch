#!/bin/bash
#SBATCH --job-name=6b.bqsrgatk
#SBATCH --output=results/slurm_logs/variants2/6b.bqsrgatk.%j.out
#SBATCH --error=results/slurm_logs/variants2/6b.bqsrgatk.%j.err
#SBATCH --time=12:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#SBATCH --ntasks-per-node 8
#SBATCH --mem=36G
#SBATCH --array=1-535

set -x
module purge
source ~/.bashrc
conda activate gatk4.2.5.0
cd $1

mkdir -p results/variants2/gatk/

##Get interval file name
intervals=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' intervals2mb.list)

##Get interval set number
intervalID=$(echo "$intervals" | cut -f2 -d_ | cut -f1 -d.)

##get bam directory
bamDir="/scratch/summit/caitlinv@colostate.edu/cawa-breed-wglc/results/bqsr"
reference="/projects/caitlinv@colostate.edu/genomes/cardellina_canadensis_pseudohap_v1.fasta"

##assign outname
mkdir -p results/variants/gatk
outname=$(printf 'results/variants2/gatk/cawa2mb.interval.bqsr.%s.vcf.gz'  $intervalID)

##call snps
gatk --java-options "-Xmx30g -XX:+UseParallelGC -XX:ParallelGCThreads=8" HaplotypeCaller \
-R $reference \
$(printf ' -I %s ' $bamDir/*recal.bam) \
-L $intervals \
-O ${outname}
