#!/bin/bash 
#SBATCH --job-name=13b.sfs.fst 
#SBATCH --output=results/slurm_logs/esu/13b.sfs.fst.%A.%a.out
#SBATCH --error=results/slurm_logs/esu/13b.sfs.fst.%A.%a.err
#SBATCH -t 01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --array=1-120
#SBATCH --mail-type=END
#SBATCH  --mail-user=caitlinv@colostate.edu


set -x

module purge
source ~/.bashrc 
cd $1
eval $(line_assign.sh $SLURM_ARRAY_TASK_ID resources/fst.combinations.1X.txt)

conda activate angsd
outDir=results/esu/ibd
name="minInd0.7.majmin4"

/projects/caitlinv@colostate.edu/software/angsd/misc/realSFS -P 8  $outDir/"$sample"."$name".saf.idx  $outDir/"$pop2"."$name".saf.idx 2> /dev/null > $outDir/$sample.$pop2.$name.sfs

###need to define number of sites
/projects/caitlinv@colostate.edu/software/angsd/misc/realSFS -P 8  $outDir/"$sample"."$name".saf.idx  $outDir/"$pop2"."$name".saf.idx   > $outDir/$sample.$pop2.$name.samesite.estim.ml

##FST for each site
/projects/caitlinv@colostate.edu/software/angsd/misc/realSFS fst index  -P 8  $outDir/"$sample"."$name".saf.idx  $outDir/"$pop2"."$name".saf.idx  -sfs $outDir/$sample.$pop2.$name.samesite.estim.ml -fstout $outDir/$sample.$pop2.$name.realsfs -whichFST 1

##FST for 50kb window, 25 kb slide
/projects/caitlinv@colostate.edu/software/angsd/misc/realSFS  fst stats2 $outDir/$sample.$pop2.$name.realsfs.fst.idx -win 50000 -step 25000 -whichFST 1 > $outDir/$sample.$pop2.$name.fst.50kbwin.txt

##FST for global
/projects/caitlinv@colostate.edu/software/angsd/misc/realSFS  fst stats $outDir/$sample.$pop2.$name.realsfs.fst.idx > $outDir/$sample.$pop2.$name.fst.global.txt
