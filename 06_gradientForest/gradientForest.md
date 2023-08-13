# CAWA gradient forest
Package versions
```
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

This is how I created the gradient forest model of genotype-environment associations and then predicted future relationships. Using the future relationships I found the euclidean distance between the current and future gene-environment associations to predict genomic offset. 

Workflow
1. Generate SAF for the adaptive loci using angsd.
2. Create a gradientForest model and test it against random chance.
3. Predict the future gene-environment relationships and measure the distance between current and future allele frequencies to find genomic offset.

## 1) Generate SAF for the adaptive loci using angsd.

Scripts: 
  1) [get.adapt.saf.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/06_gradientForest/scripts/get.adapt.saf.sbatch)
  2) [maf.combine.R](https://github.com/mcaitlinv/cawa-breeding/blob/main/06_gradientForest/scripts/maf.combine.R)

To use gradient forest for gene-environment relationships we need to get a site allele frequency for each SNP in the adaptive loci set for each sampled site. I used angsd to calculate this from genotype likelihoods. Once I had an SAF for each site, I combined all sites into a single dataframe. From this I subset the data to only SNPs with no missing data. 

## 2) Create a gradientForest model and test it against random chance.

Scripts: 
  1) [gradientForestRandomization.R](https://github.com/mcaitlinv/cawa-breeding/blob/main/06_gradientForest/scripts/gradientForestRandomization.R)

Here I used the environmental variables as predictors with the allele frequencies as the response. This generates a turnover function for each SNP that corresponds to the environmental predictors. Once we have this, to check the model is better than random chance, I randomized the predictors and created 100 random models. Those models, when compared to our initial model, were worse at predicting the associations based on R-squared. 

## 3) Predict the future gene-environment relationships and measure the distance between current and future allele frequencies to find genomic offset.

Scripts: 
 1) [gradForAdaptive.Rmd]((https://htmlpreview.github.io/?https://github.com/mcaitlinv/cawa-breeding/blob/main/06_gradientForest/scripts/gradForAdaptive.html)
  
After creating a model that could predict SNP turnover across the current breeding range, I used the model to predict allele frequency- assuming gene-environment relationships stayed the same- in 2060-2080. To do this, I used the predicted BIOCLIM values from 100,000 random points across the breeding range in 2060-2080 according to shared socio-economic pathways(SSP) 126 and 585. The values in these pathways represent the 'best' and 'worst' case scenarios of predicted climate change from 2060-2080. Using the model with BIOCLIM values as predictors, I found the predicted allele frequencies needed to maintain the same gene-environment relationship given environmental change. Once I had the predicted 'future' allele frequency set, I calculated the euclidean distance between current allele frequencies and the 'future' set to find genomic offset. The larger the value the more predicted change is necessary to maintain gene-environment relationships. I then visualized this across the breeding range and estimated the median genomic offset across the breeding range and within adaptive units.
