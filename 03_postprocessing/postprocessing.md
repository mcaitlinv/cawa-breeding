# Post processing
Package versions
```
bcftools 1.15.1
GATK 4.2.5.0
picard  2.26.11
samtools  1.11
bedtools 2.30.0
beagle 4.1
R version 4.1.2
R packages:
  srsStuff 0.1.0
  tidyverse 1.3.1
  RColorBrewer 1.1-3
```

These steps go from the vcf to the post-processed vcf. I had to do some post-processing before using these variants in future analyses due to sequencing bias.

Workflow
1. Initial check by PCA for any clustering due to sex, coverage, etc. 
2. Find and remove variants based on sequencing platform
3. Re-check PCA
4. Impute missing data

## 1) Initial check by PCA for any clustering due to sex, coverage, etc.

Scripts: 
  1) [pca.sequencer.rmd](https://htmlpreview.github.io/?https://github.com/mcaitlinv/cawa-breeding/blob/main/03_postprocessing/scripts/pca.sequencer.html)
  
I used single read sampling to even out coverage across bams, since low-coverage samples can have read-depth bias that looks like population structure. I used a mix of rmd scripts and data manipulations using R on the HPCC using a script called runR.sbatch to submit R scripts with extra memory. 

## 2) Find and remove variants based on sequencing platform

Scripts: 
  1) [get.FSTbams.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/03_postprocessing/scripts/get.FSTbams.sbatch)
  2) [get.conda.fstgatk.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/03_postprocessing/scripts/get.conda.fstgatk.sbatch)
  3) [get.mergevcfs.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/03_postprocessing/scripts/get.fstvcf.sbatch)
  4) [get.fstFilter.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/03_postprocessing/scripts/get.fstFilter.sbatch)
  5) [get.fstvcf.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/03_postprocessing/scripts/get.fstvcf.sbatch)

After finding platform effects, I decided to split samples that were run on both platforms into platform-based bams, then find FST values within each sample. The logic being that the same individuals would have highest FST values based on platform differences. 

  - I split bams using `samtools split`, then processed bams using GATK with 2Mb section to find variants. 
  - I merged the sections back together in a single vcf.
  - I then found FST and calculated the FST percentiles across the groups. I took the 80th, 85th, 90th, 95th, and 99th percentiles and pulled out the variants with higher FST values based on `vcftools` FST pairwise calculation. 
  - Then I filtered out variants in the full dataset that corresponded with the FST percentiles and generated the allele depth file to check these variant sets via PCA.

## 3) Re-check PCA

Scripts: 
  1) [pca.sequencer.rmd](https://htmlpreview.github.io/?https://github.com/mcaitlinv/cawa-breeding/blob/main/03_postprocessing/scripts/pca.sequencer.html)
  
Using the filtered variant sets, I used PCA to check when platform associated clustering on PCA was not the primary separation. I decided the 85th percentile group had the least platform effects while still maintaining as many variants as possible. I used this variant set as my filtered, primary variant set for the rest of the analysis. 

## 4) Impute missing data

Scripts: 
  1) [8a.conda.beagle.impute.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/03_postprocessing/scripts/8a.conda.beagle.impute.sbatch)
  2) [8b.conda.beagle.impute.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/03_postprocessing/scripts/8b.conda.beagle.impute.sbatch)
  3) [8c.conda.beagle.impute.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/03_postprocessing/scripts/8c.conda.beagle.impute.sbatch)

Because some of the downstream analyses do not deal with missing data very well, I imputed missing data using `beagle`. Initial filtering filtered out any variant with more than 20% missing across individuals. Imputation was done using genotype likelihoods. 
