#!/bin/bash
#set a job name
#SBATCH --job-name=getFSTGATK
#SBATCH --output=results/slurm_logs/postprocess/getFSTGATK.%j.out
#SBATCH --error=results/slurm_logs/postprocess/getFSTGATK.%j.err
################
#SBATCH --time=4:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#################
# Note: 4.84G/core or task
#################
#SBATCH --ntasks-per-node 8
#SBATCH --mem=36G
#SBATCH --array=1-535
#################
set -x

module purge
source ~/.bashrc
conda activate gatk4.2.5.0
cd $1

##Get interval file name
intervals=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' intervals2mb.list)
##Get interval set number
intervalID=$(echo "$intervals" | cut -f2 -d_ | cut -f1 -d.)

##get bam directory
bamDir="/scratch/summit/caitlinv@colostate.edu/cawa-breed-wglc/results/fstbams"
reference="/projects/caitlinv@colostate.edu/genomes/cardellina_canadensis_pseudohap_v1.fasta"

##assign outname
mkdir -p results/fstbams/vcfs
outname=$(printf 'results/fstbams/vcfs/cawa2mb.interval.%s.vcf.gz'  $intervalID)

gatk --java-options "-Xmx30g -XX:+UseParallelGC -XX:ParallelGCThreads=8" HaplotypeCaller \
-R $reference \
$(printf ' -I %s ' $bamDir/*newRG.bam) \
-L $intervals \
-O ${outname}
