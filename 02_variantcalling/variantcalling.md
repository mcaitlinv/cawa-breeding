# Variant calling

These steps go from pre-processed bams through variant calling and filtering to get the dataset of SNPs to analyze. 

Workflow
4. Initial variant calling with GATK and mpileup: 
5. Merge/filter and intersect variant sets:
6. Base quality score recalibration:
7. Final variant calling:

Each step of the workflow has 2 associated scripts, eg scripts 4a/4b are in the first step.

## 4) Initial variant calling with GATK and mpileup

To get a list of 'known' variants to perform BQSR with, I used GATK and mpileup to call SNPs with the pre-processed bams. Both of these steps can take a loooong time if you try to use the full scaffold set at once. Instead, I split my scaffolds into regions of ~2Mb each and called variants for each 2Mb region. To get the 2Mb regions, I used the script `get.2MB.scaffolds.sbath`, which uses the bam header to calculate scaffold size, then splits scaffolds that are larger than 1Mb into ~1Mb chunks, then adds scaffolds together until it reaches ~2Mb. This results in both lots of scaffolds that are tiny in a file and some files with just 2 1Mb chunks of a single scaffold. For example:

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
  
  -Then I used `GATK HaplotypeCaller` to get a file of called snps. 



