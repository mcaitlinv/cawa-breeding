#!/bin/bash
#SBATCH --job-name=13a.saf
#SBATCH --output=results/slurm_logs/esu/13a.saf.%A.%a.out
#SBATCH --error=results/slurm_logs/esu/13a.saf.%A.%a.err
#SBATCH --partition=shas
#SBATCH --qos=normal
#SBATCH -t 06:00:00
#SBATCH --nodes=1
#SBATCH --array=1-16
#SBATCH --ntasks=6
#SBATCH --mail-type=END
#SBATCH  --mail-user=caitlinv@colostate.edu

set -x
module purge
source ~/.bashrc
cd $1
eval $(line_assign.sh $SLURM_ARRAY_TASK_ID resources/populations.bamlist.list) 
outfile=`basename $sample .bamlist`.minInd0.7
num=$(wc -l $sample | awk '{printf "%.0f", $1*0.7}')
outDir=results/esu/ibd 

conda activate angsd

##Make sites to use for pairwise fst comparisons
##Uses previously filtered sites from vcf minus any platform effects
##run once
# sort -k 1 -n -b -f -t $'\t' resources/sites.filtRda.txt > resources/sites.sort.filtRda.txt
# angsd sites index resources/sites.sort.filtRda.txt

# ##get unique scaffolds
# cut -d' ' -f1 resources/sites.sort.filtRda.txt | sort -u > resources/angsd.scaffs.txt

#get saf for fst 
/projects/caitlinv@colostate.edu/software/angsd/angsd -bam $sample -P 6 \
	-ref /projects/caitlinv@colostate.edu/genomes/cardellina_canadensis_pseudohap_v1.fasta \
	-anc /projects/caitlinv@colostate.edu/genomes/cardellina_canadensis_pseudohap_v1.fasta \
	-rf resources/angsd.scaffs.txt \
	-dosaf 1 \
	-uniqueOnly 1 -remove_bads 1 -baq 1 -C 50 \
	-doGlf 2 \
	-minInd $num -doCounts 1 \
	-GL 1 -doMajorMinor 4  \
	-doGeno 32 -doPost 1 \
	-sites resources/sites.sort.filtRda.txt -out $outDir/"$outfile".majmin4
