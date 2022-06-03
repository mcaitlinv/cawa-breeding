# Variant calling

These steps go from pre-processed bams through variant calling and filtering to get the dataset of SNPs to analyze. 

Workflow
1. Initial variant calling with GATK and mpileup: 
2. Merge/filter and intersect variant sets:
3. Base quality score recalibration:
4. Final variant calling:

Each step of the workflow has 2 associated scripts, eg scripts 4a/4b are in the first step.

## 1) Initial variant calling with GATK and mpileup

To get a list of 'known' variants to perform BQSR with, I used GATK and mpileup to call SNPs with the pre-processed bams. Both of these steps can take a loooong time if you try to use the full scaffold set at once. Instead, I split my scaffolds into regions of ~2Mb each and called variants for each 2Mb region. To get the 2Mb regions, I used the script `get.2MB.scaffolds.sbatch`, which uses the bam header to calculate scaffold size, then splits scaffolds that are larger than 1Mb into ~1Mb chunks, then adds scaffolds together until it reaches ~2Mb. This results in both lots of scaffolds that are tiny in a file and some files with just 2 1Mb chunks of a single scaffold. For example:

``` 
cat intervals/interval2mb_003.list
scaffold4:4000001-5000000
scaffold4:5000001-6000000

cat intervals/interval2mb_473.list
scaffold29536:1000001-1735654
scaffold29537:1-1238
scaffold29538:1-1585
scaffold29539:1-1327
scaffold29540:1-9581
scaffold29541:1-8376
scaffold29542:1-26494
scaffold29543:1-2448
scaffold29544:1-5374
scaffold29545:1-31126
scaffold29546:1-1000000
scaffold29546:1000001-2000000

```
  
  - To get snps I first used `bcftools mpileup` to get genotype likelihoods, then the `bcftools call` function to produce a file with called snps.
  
  - Then I used `GATK HaplotypeCaller` to get a file of called snps. 
  
Mpileup variant calling:

``` bash

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
```

GATK variant calling:

``` bash

source ~/.bashrc
conda activate gatk4.2.5.0
cd $1

##Get interval file name
intervals=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' intervals2mb.list)
##Get interval set number
intervalID=$(echo "$intervals" | cut -f2 -d_ | cut -f1 -d.)

##get bam directory
bamDir="/scratch/summit/caitlinv@colostate.edu/cawa-breed-wglc/results/bams"
reference="/projects/caitlinv@colostate.edu/genomes/cardellina_canadensis_pseudohap_v1.fasta"

##assign outname
mkdir -p results/variants/gatk
outname=$(printf 'results/variants/gatk/cawa2mb.interval.%s.vcf.gz'  $intervalID)

gatk --java-options "-Xmx30g -XX:+UseParallelGC -XX:ParallelGCThreads=8" HaplotypeCaller \
-R $reference \
$(printf ' -I %s ' $bamDir/*mergeMkDup.bam) \
-L $intervals \
-O ${outname}
```

## 2) Merge and filter/intersect variant sets:

Once I have my initial vcfs generated with `bcftools` and `GATK`, I merge the scaffolds together into a single vcf then filter and intersect the variants called with samtools and GATK. This is to produce a 'high quality' variant set per the BQSR best practices for non-model references, see ['No excuses'](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-). I do with with slightly more stringent filters than I would normally since the goal is to have the best possible variant set for BQSR.

  - To merge the scaffold vcfs into a single vcf, I use `bcftools index` to index, then use `bcftools concat` to get the vcfs into a single vcf. Note that I used the -a option in `bcftools concat` to allow positionsto be out of order when concatenation. They shouldn't be due to how I split my scaffolds, but just in case I allow them to be out of order and then use `bcftools sort` to sort the final vcf.
  
  - Then I used `bcftools view` to hard filter each set of variants. I kept only biallelic sites and set a minor allele frequncy of less than 5%. I used a missing filter of less than 5% and a quality score of above 50. This is a bit more stringent than I would normally filter to get the highest confidence variants. I then used `bcftools isec` to get the intersection of the filtered sets of SNPs. This is then used for the 'known variant set' for BQSR.
  
Merging vcfs:

``` bash 
source ~/.bashrc
cd $1

conda activate bcftools

##Make temp dir, bcftools was not finding default tmp
mkdir -p tmp

# list file with the 535 2mb intervals
vcfs=intervals2mb.vcfs.list
bcfs=intervals2mb.bcfs.list

##out names for merged files
outname=cawa.unfiltered

##index bcfs
while read a; do bcftools index $a; done < intervals2mb.bcfs.list
##bcftools concat for bcfs
bcftools concat -a  --threads 8 -f ${bcfs} -Oz -o results/variants/"$outname".bcf.gz
bcftools sort -m 30G -T tmp -Oz -o results/variants/"$outname".sort.bcf.gz  results/variants/"$outname".bcf.gz

##No need to index, gatk already made it
##bcftools concat for bcfs
bcftools concat -a  --threads 8 -f ${vcfs} -Oz -o results/variants/"$outname".vcf.gz
bcftools sort -m 30G -T tmp -Oz -o results/variants/"$outname".sort.vcf.gz  results/variants/"$outname".vcf.gz
```
Filtering and intersecting vcfs:

``` bash 
source ~/.bashrc
cd $1

conda activate bcftools

##Make temp dir, bcftools was not finding default tmp
mkdir -p tmp

##the 2 vcfs to filter then intersect
vcf=temp.sort.vcf.gz
bcf=temp.sort.bcf.gz
##out directory
outDir="results/variants"

##index the vcfs 
bcftools index $outDir/$vcf
bcftools index $outDir/$bcf
##Get unfiltered stats
bcftools stats $outDir/$vcf > $outDir/$vcf.stats
bcftools stats $outDir/$bcf > $outDir/$bcf.stats

##filter vcfs
outname=cawa.snps.miss0.1.qual30
bcftools view --threads 8 -m 2 -M 2 --min-af 0.05 --max-af 0.95 -i 'F_MISSING < 0.1 && QUAL>30' -Oz -o $outDir/$outname.vcf.gz $outDir/$vcf
bcftools view --threads 8 -m 2 -M 2 --min-af 0.05 --max-af 0.95 -i 'F_MISSING < 0.1 && QUAL>30' -Oz -o $outDir/$outname.bcf.gz $outDir/$bcf
##index filtered
bcftools index $outDir/"$outname".vcf.gz
bcftools index $outDir/"$outname".bcf.gz
##stats on filtered
bcftools stats $outDir/"$outname".vcf.gz > $outDir/"$outname".vcf.stats
bcftools stats $outDir/"$outname".bcf.gz > $outDir/"$outname".bcf.stats

##intersect variant sets, get stats on intersection
mkdir -p $outDir/isec
bcftools isec $outDir/"$outname".vcf.gz $outDir/"$outname".bcf.gz -p $outDir/isec -Oz
bcftools index $outDir/isec/*
bcftools stats $outDir/isec/0002.vcf.gz > $outDir/isec/0002.vcf.stats
```


