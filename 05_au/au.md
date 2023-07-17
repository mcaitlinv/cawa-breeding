# CAWA AU
Package versions
```
bcftools 1.15.1
admixture 1.3.0
plink 1.9
angsd 
R version 4.1.2
R packages:
  rgdal
  extendedForest
  gradientForest
  gstat
  RColorBrewer
  rasterVis
  tidyverse
  readr
  psych
  sf
  ggpub
  raster
  cluster
  ggplot2
  readxl
  smoothr
  lwgeom
  rgeos 
```

These steps go from the post-filtering vcf to putatively adaptive loci. I used the putatively adaptive loci to look for adaptive structure and then visualize the structure across the CAWA breeding range. 

Workflow
1. Use gradientForest to identify most important environmental variables.
2. Use redundancy analysis and the environmetal variables to find putatively adaptive loci.
3. Use latent factor mixed models and the environmetal variables to find putatively adaptive loci.
4. Identify adaptive loci structure using ADMIXTURE.
5. Visual the adaptive loci in a spatially explicit manner.

## 1) Use gradientForest to identify most important environmental variables.

Scripts: 
  1) [envVar.rmd](https://github.com/mcaitlinv/cawa-breeding/blob/main/05_au/scripts/envVar.rmd)
  2) [gradFor_multiRun.R](https://github.com/mcaitlinv/cawa-breeding/blob/main/05_au/scripts/gradFor_multiRun.R)

I used BIOCLIM variables, tree cover, elevation, and a human influence index as environmental variable set. I extracted the variables for each site and individual coordinate. Using BIOCLIM predictions for future climate in 2060-2080, I also extracted variables for future cliamte and generated 100,000 random points across the CAWA breeding range. Using only the current day environmental variables, I ran gradientForest with 5 subsets of 50,000 SNPs from the dataset. This was to identify what variables were important to genetic diversity throughout the genome. I identified 4 uncorrelated variables- BIO

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
