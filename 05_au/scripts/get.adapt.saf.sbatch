#!/bin/bash
#SBATCH --job-name=adapt.saf
#SBATCH --output=results/slurm_logs/au/adapt.saf.%A.%a.out
#SBATCH --error=results/slurm_logs/au/adapt.saf.%A.%a.err
#SBATCH --partition=shas
#SBATCH --qos=normal
#SBATCH -t 12:00:00
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
outfile=`basename $sample .bamlist`.minInd0.5
num=$(wc -l $sample | awk '{printf "%.0f", $1*0.5}')
outDir=results/au/gradientForAdaptOnly

conda activate angsd

##Make sites to use for pairwise fst comparisons
##Uses previously filtered sites from vcf minus any platform effects
##run once
# zgrep -v '#' results/au/cawa.adapt.all.vcf.gz | awk '{print $1, $2}' > resources/adapt.sites.txt
# sort -k 1 -n -b -f -t $'\t' resources/adapt.sites.txt > resources/adapt.sort.sites.txt
# angsd sites index resources/adapt.sort.sites.txt

# # ##get unique scaffolds
# cut -d' ' -f1 resources/adapt.sort.sites.txt | sort -u > resources/angsd.adapt.scaffs.txt

#get saf for fst 
/projects/caitlinv@colostate.edu/software/angsd/angsd -bam $sample -P 6 \
	-ref /projects/caitlinv@colostate.edu/genomes/cardellina_canadensis_pseudohap_v1.fasta \
	-anc /projects/caitlinv@colostate.edu/genomes/cardellina_canadensis_pseudohap_v1.fasta \
	-rf resources/angsd.adapt.scaffs.txt \
	-dosaf 1 \
	-uniqueOnly 1 -remove_bads 1 -baq 1 -C 50 \
	-doGlf 2 \
	-minInd $num -doCounts 1 \
	-GL 1 -doMajorMinor 4 -doMaf 1 \
	-doGeno 32 -doPost 1 \
	-sites resources/adapt.sort.sites.txt -out $outDir/"$outfile".majmin4.adapt
