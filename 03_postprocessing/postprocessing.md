# Post processing
Package versions
```
bcftools 1.15.1
GATK 4.2.5.0
picard  2.26.11
samtools  1.11
bedtools 2.30.0
R version 4.1.2
R packages:
  srsStuff 0.1.0
  tidyverse 1.3.1
  RColorBrewer 1.1-3
```

These steps go from the vcf to the post-processed vcf. I had to do some post-processing before using these variants in future analyses due to sequencing bias.

Workflow
1. Initial check by PCA for any clustering due to sex, coverage, etc. 
2. Find variants based on sequencing platform
3. Remove variants based on sequencing platform
4. Re-check PCA

## 1) Initial variant calling with GATK and mpileup

Scripts: 
  1) [4a.conda.snpcall.mpileup.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/02_variantcalling/scripts/4a.conda.snpcall.mpileup.sbatch)
  2) [4b.conda.snpcall.gatk.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/02_variantcalling/scripts/4b.conda.snpcall.gatk.sbatch)

