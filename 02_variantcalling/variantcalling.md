# Variant calling
Package versions
```
bcftools 1.15.1
GATK 4.2.5.0
```

These steps go from pre-processed bams through variant calling and filtering to get the dataset of SNPs to analyze. 

Workflow
1. Initial variant calling with GATK and mpileup: 
2. Merge/filter and intersect variant sets:
3. Base quality score recalibration:
4. Variant filtering:

Each step of the workflow has 2 associated scripts, eg scripts 4a/4b are in the first step.

## 1) Initial variant calling with GATK and mpileup

Scripts: 
  1) [4a.conda.snpcall.mpileup.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/02_variantcalling/scripts/4a.conda.snpcall.mpileup.sbatch)
  2) [4b.conda.snpcall.gatk.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/02_variantcalling/scripts/4b.conda.snpcall.gatk.sbatch)

To get a list of 'known' variants to perform BQSR with, I used GATK and mpileup to call SNPs with the pre-processed bams. Both of these steps can take a loooong time if you try to use the full scaffold set at once. 
Instead, I split my scaffolds into regions of ~2Mb each and called variants for each 2Mb region. To get the 2Mb regions, I used the script `get.2MB.scaffolds.sbatch`, which uses the bam header to calculate scaffold size, then splits scaffolds that are larger than 1Mb into ~1Mb chunks, then adds scaffolds together until it reaches ~2Mb. This results in both lots of scaffolds that are tiny in a file and some files with just 2 1Mb chunks of a single scaffold. For example:

```bash
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
  
  - To get snps I first used `bcftools mpileup` to get genotype likelihoods, then the `bcftools call` function to produce a file with called snps for each 2Mb region.
  
  - Then I used `GATK HaplotypeCaller` to get a multisample vcf of called snps for each 2Mb region.
  
## 2) Merge and filter/intersect variant sets:

Scripts: 
  1) [5a.conda.mergeVcfs.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/02_variantcalling/scripts/5a.conda.mergeVcfs.sbatch)
  2) [5b.conda.filter.intersectVcfs.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/02_variantcalling/scripts/5b.conda.filter.intersectVcfs.sbatch)

Once I have my initial vcfs generated with `bcftools` and `GATK`, I merge the scaffolds together into a single vcf then filter and intersect the variants called with bcftools and GATK. This is to produce a 'high quality' variant set per the BQSR best practices for non-model references, see ['No excuses'](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-). I do with with slightly more stringent filters on missingness than I would normally since the goal is to have the best possible variant set for BQSR.

  - To merge the scaffold vcfs into a single vcf, I use `bcftools index` to index, then use `bcftools concat` to get the vcfs into a single vcf. Note that I used the -a option in `bcftools concat` to allow positions to be out of order when concatenating. They shouldn't be due to how I split my scaffolds, but just in case I allow them to be out of order and then use `bcftools sort` to sort the final vcf.
  
  - Then I used `bcftools view` to hard filter each set of variants. I kept only bi-allelic sites and set a minor allele frequency of less than 5%. I used a missing filter of less than 10% and a quality score of above 30. I then used `bcftools isec` to get the intersection of the filtered sets of SNPs. This is then used for the 'known variant set' for BQSR.

## 3) Base quality score recalibration:

Scripts: 
  1) [6a.conda.bqsr.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/02_variantcalling/scripts/6a.conda.bqsr.sbatch)
  2) [6b.conda.recalgatk.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/02_variantcalling/scripts/6b.conda.recalgatk.sbatch)

After getting a 'high quality' variant set, I used base quality score recalibration to make a recalibrated bam for each sample. I used the new recalibrated bams to call a new variant set. 

  - With the 'high quality' variant set, I used `GATK BaseRecalibrator` to get recalibration tables for each sample. Then using the reaclibration tables, I used `GATK ApplyBQSR` to get recalibrated bams for each sample. I generated recalibration tables from the recalibrated bams using `GATK BaseRecalibrator` and then used `GATK AnalyzeCovariates` to see the differences in base quality before and after recalibration. 
  
  - Then I used `GATK HaplotypeCaller` to get a multi-sample vcf for each of the regions of 2Mb lengths of scaffold.

## 4) Variant filtering:

Scripts: 
  1) [7a.conda.mergeVcfs.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/02_variantcalling/scripts/7a.conda.mergeVcfs.sbatch)
  2) [7b.conda.filtervcfs.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/02_variantcalling/scripts/7b.conda.filtervcfs.sbatch)

Once I had the vcf files for all of the 2Mb regions, I merged them into a single large vcf. Then I filtered the vcf to get a final set of SNPs.

  -  To merge the scaffold vcfs into a single vcf, I used `bcftools index` to index, then use `bcftools concat` to get the vcfs into a single vcf. Note that I used the -a option in `bcftools concat` to allow positions to be out of order when concatenating. I used `bcftools sort` to sort the final vcf.
  
  -  To get the final variant set, I used `bcftools view` to hard filter the vcf. I kept only bi-allelic sites and set a minor allele frequency of less than 5%. I used a missing filter of less than 20% and a quality score of above 30.


