# CAWA ESU
Package versions
```
admixture 1.3.0
plink 1.9
angsd 
eems
R version 4.1.2
R packages:
  readr
  tidyverse
  ggplot2
  gridExtra
  pophelper
  gdata
  readxl
  srsStuff
  RColorBrewer
```

These steps go from the post-filtering vcf to population structure analyses. I delineated ESUs using ADMIXTURE, PCA, and EEMS. 

Workflow
1. Use ADMIXTURE to find suggested K value for the entire population.
2. Identify clusters on PCA of whole genome for population structure.
3. Check for isolation by distance using pairwise FST between sites.
4. Use EEMS to check for areas of lower than expected gene flow to identify potential barriers.

## 1) Use ADMIXTURE to find suggested K value for the entire population.

Scripts: 
  1) [9.conda.admixture.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/04_esu/scripts/9.conda.admixture.sbatch)

Data visualization:
 1) [admixture.rmd](https://htmlpreview.github.io/?https://github.com/mcaitlinv/cawa-breeding/blob/main/04_esu/scripts/admixture.html)
 
I filtered out variants in linkage disequilibrium  using `plink` and ran 5 random seeds of `ADMIXTURE` across K's 1-6. I used the estimates of CV error to find the best K.

## 2) Identify clusters on PCA of whole genome for population structure.

Scripts: 
  1) [SRS_PCA.rmd](https://htmlpreview.github.io/?https://github.com/mcaitlinv/cawa-breeding/blob/main/04_esu/scripts/SRS_PCA.html)
  
I used single read sampling PCA to check for population structure. Single read sampling was appropriate because using low coverage individuals bias can be introduced with varying read depths that looks like population structure. To avoid the spurious clustering, single read sampling evens out the coverage across individuals. 


## 3) Check for isolation by distance using pairwise FST between sites.

Scripts: 
  1) [10a.saf.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/04_esu/scripts/10a.saf.sbatch)
  2) [10b.conda.sfs.fst.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/04_esu/scripts/10b.conda.sfs.fst.sbatch)
  
Data visualization:
 1) [IBD.rmd](https://htmlpreview.github.io/?https://github.com/mcaitlinv/cawa-breeding/blob/main/04_esu/scripts/IBD.html)
  
I used `ANGSD` to get the genotype likelihood site allele frequency and then used the SAF to generate pairwise FST between the sampled sites. Using the pairwise FST between sampled sites and site location, I used mantel tests to check for isolation by distance. I also extracted environmental data at each site to use for a partial mantel test that I used to check for isolation by environment. 

## 4) Use EEMS to check for areas of lower than expected gene flow to identify potential barriers.

Scripts: 
  1) [11.conda.eems.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/04_esu/scripts/11.conda.eems.sbatch)

Data visualization:
 1) [eems.rmd](https://htmlpreview.github.io/?https://github.com/mcaitlinv/cawa-breeding/blob/main/04_esu/scripts/eems.html)
 
I used `plink` to manipulate the vcf file into EEMS format, then used `bed2diffs` that is part of `EEMS` to make a genetic dissimilarity matrix. Using the matrix and a file that had coordinates of the outer edge of the CAWA breeding range, I used EEMS to generate a estimate of the gene 'migration' surface across the breeding range. Using that data, I then rasterized and plotted the values to the CAWA breeding range.
