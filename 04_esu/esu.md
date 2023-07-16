# CAWA ESU
Package versions
```
admixture 1.3.0
plink 1.9

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
 1) [admixture.rmd](https://htmlpreview.github.io/ https://github.com/mcaitlinv/cawa-breeding/blob/main/04_esu/scripts/admixture.html)
 
I filtered out variants in linkage disequilibrium  using `plink` and ran 5 random seeds of `ADMIXTURE` across K's 1-6. I used the estimates of CV error to find the best K.

## 2) Identify clusters on PCA of whole genome for population structure.

Scripts: 
  1) [SRS_PCA.rmd](https://htmlpreview.github.io/ https://github.com/mcaitlinv/cawa-breeding/blob/main/04_esu/scripts/SRS_PCA.html)
  
I used single read sampling PCA to check for population structure. Single read sampling was appropriate because using low coverage individuals bias can be introduced with varying read depths that looks like population structure. To avoid the spurious clustering, single read sampling evens out the coverage across individuals. 


## 3) Re-check PCA

Scripts: 
  1) [pca.sequencer.rmd](https://github.com/mcaitlinv/cawa-breeding/blob/main/03_postprocessing/scripts/pca.sequencer.html)
  
Using the filtered variant sets, I used PCA to check when platform associated clustering on PCA was not the primary separation. I decided the 85th percentile group had the least platform effects while still maintining as many variants as possible. I used this variant set as my filtered, primary variant set for the rest of the analysis. 
