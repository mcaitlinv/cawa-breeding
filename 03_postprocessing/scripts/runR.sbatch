#!/bin/bash
#SBATCH --job-name=runR
#SBATCH --output=results/slurm_logs/au/adaptonly.gradFor.%j.out
#SBATCH --error=results/slurm_logs/au/adaptonly.gradFor.%j.err
#SBATCH --partition=shas
#SBATCH --qos=normal
#SBATCH -t 12:00:00
#SBATCH --mem=24G
#SBATCH --ntasks 4
#SBATCH --mail-type=END
#SBATCH  --mail-user=caitlinv@colostate.edu

set -x

module purge
eval "$('/projects/caitlinv@colostate.edu/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
cd $1
conda activate r_env
name=`basename "$2"`
outDir=results/au/gradientForAdaptOnly

R CMD BATCH  $2 $outDir/"$name".221011.Rout
##for gradient forest
#R CMD BATCH '--args '$3'' $2 results/au/gradientFor/"$name"."$3".random.220703.Rout
##for conStruct
#R CMD BATCH '--args '$SLURM_ARRAY_TASK_ID'' $2 $outDir/"$name".220831."$SLURM_ARRAY_TASK_ID".Rout

#cat $outDir/cawa.adaptOnly.16ClimGroup.minInd0.5.majmin4.NA.txt| cut -f 1-18 | sed 1d | awk 'BEGIN{FS="\t";OFS="\t"};{count=0;for(i=3; i<19; i++) {if($i!="NA") {count++}};if (count>=15){print $0}}' > $outDir/cawa.adaptOnly.16ClimGroup.minInd0.5.majmin4.scaf.noNA.txt



