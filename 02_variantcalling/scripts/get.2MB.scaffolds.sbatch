#!/bin/bash 
#set a job name  
#SBATCH --job-name=get2MBscaffolds
#SBATCH --output=results/slurm_logs/variants/get2MBscaffolds.%j.out
#SBATCH --error=results/slurm_logs/variants/get2MBscaffolds.%j.err
#################
#SBATCH -t 00:05:00
#SBATCH --partition=shas
#SBATCH --nodes=1
#SBATCH --mem=4G
#################
#SBATCH --mail-type=END
#SBATCH  --mail-user=caitlinv@colostate.edu
#################

set -x

module purge
source ~/.bashrc
cd $1
conda activate bwasam

##make directory to store things
mkdir -p intervals

##Get scaffold names and lengths from bam header
samtools view  -H results/bams/14N1963.mergeMkDup.bam | sed 's/SN://g; s/LN://g' | awk '/^@SQ/' | awk '{print $2, $3}' > intervals/scaffold.lengths

##Split scaffolds into 1MB or less segments and get length of segment in column 2
cat intervals/scaffold.lengths | awk '{tot=0; while(tot<$2) {start=tot+1; tot+=1e6; if(tot>$2) tot=$2; printf("%s\t%d\t%d\t%d\n", $1,start,tot, (tot-start+1));}}' > intervals/intervals1mb.txt

##Move to intervals directory since this awk command gets testy
cd intervals

##Add the scaffolds column 2 together until it reaches 2MB, then print the scaffolds to a file
awk 'BEGIN {id=1; tot=0} {if(tot<2e6) {tot+=$4; k=sprintf("%03d",id); print $1,$2,$3 >> "interval2mb_" k ".txt"} else {id+=1; tot=$4; k=sprintf("%03d", id); print $1,$2,$3 >> "interval2mb_" k ".txt";}}' intervals1mb.txt

##Now get it in the GATK intervals format
for f in intervals/interval2mb*.txt;
do 
awk '{printf("%s:%d-%d\n", $1,$2,$3)}' $f > `basename $f .txt`.list
done

##Get a list of interval file names
cd $1
ls intervals/interval2mb*.txt > intervals2mb.txt
ls intervals/interval2mb*.list > intervals2mb.list
