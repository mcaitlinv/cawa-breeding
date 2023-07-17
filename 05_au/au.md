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
  vegan
  lfmm
  qvalue
  adegenet
  data.table
  reshape2
  codep
  adespatial
  adegraphics
```

These steps go from the post-filtering vcf to putatively adaptive loci. I used the putatively adaptive loci to look for adaptive structure and then visualize the structure across the CAWA breeding range. 

Workflow
1. Use gradientForest to identify most important environmental variables.
2. Use redundancy analysis and latent factor mixed models to find putative adaptive loci.
3. Identify adaptive loci structure using ADMIXTURE.
4. Visualize the adaptive loci in a spatially explicit manner.

## 1) Use gradientForest to identify most important environmental variables.

Scripts: 
  1) [envVar.rmd](https://github.com/mcaitlinv/cawa-breeding/blob/main/05_au/scripts/envVar.rmd)
  2) [gradFor_multiRun.R](https://github.com/mcaitlinv/cawa-breeding/blob/main/05_au/scripts/gradFor_multiRun.R)

I used BIOCLIM variables, tree cover, elevation, and a human influence index as environmental variable set. I extracted the variables for each site and individual coordinate. Using BIOCLIM predictions for future climate in 2060-2080, I also extracted variables for future cliamte and generated 100,000 random points across the CAWA breeding range. Using only the current day environmental variables, I ran gradientForest with 5 subsets of 50,000 SNPs from the dataset. This was to identify what variables were important to genetic diversity throughout the genome. I identified 4 uncorrelated variables-
    Rank1 Bio15: precipitation seasonality
    Rank5 Bio13: precipitation of the wettest month
    Rank6 Bio10: mean temperature of the warmest quarter
    Rank13 Tree: tree cover %

## 2) Use redundancy analysis and latent factor mixed models to find putative adaptive loci.

Scripts: 
  1) [lfmm.221123.rmd](https://htmlpreview.github.io/?https://github.com/mcaitlinv/cawa-breeding/blob/main/05_au/scripts/lfmm.221123.html)
  2) [rda.all.221114.rmd](https://htmlpreview.github.io/?https://github.com/mcaitlinv/cawa-breeding/blob/main/05_au/scripts/rda.all.221114.html)
  
To identify putative adaptive loci, I used both LFMM and RDA to identify potential variants. LFMM is a univariate approach that controls for population structure with latent factors, while RDA is a multivariate constrained ordination approach that performs better at finding many loci of small effect (Forester et al., 2018). Both approaches find SNPs associated with differences in environment at each individual’s given latitude and longitude. For LFMM I estimated neutral population structure to be K=3 to account for any structure we missed in the dataset, while for RDA I generated Moran's Eigen Vectors as proxy for geographic distance. 

- Note that the RDA script includes the first time I ran the RDA which did find platform effects as a major axis. I removed those effects from the dataset before doing any analysis. 

## 3) Identify adaptive loci structure using ADMIXTURE.

Scripts: 
  1) [get.adaptfilter.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/05_au/scripts/get.adaptfilter.sbatch)
  2) [get.adapt.admix.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/05_auu/scripts/get.adapt.admix.sbatch)
  3) [get.adapt.highgrading.bias.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/05_auu/scripts/get.adapt.highgrading.bias.sbatch)
  
Data visualization:
 1) [admixture.rmd](https://htmlpreview.github.io/?https://github.com/mcaitlinv/cawa-breeding/blob/main/04_esu/scripts/admixture.html)
  
Using the SNPs found with LFMM and RDA, I filtered SNPs from the dataset if they were not putatively adaptive using `bcftools`. Then I used `plink` to format the data and ran 5 random seeds of `ADMIXTURE` across K's 1-6. I used the estimates of CV error to find the best K. Additionally since there was some concern that by choosing variants and then predicting K with the entire set of samples would bias the results due to highgrading bias, I randomized the samples used for a training and test dataset. 

## 4) Visualize the adaptive loci in a spatially explicit manner.

Scripts: 
  1) [genoscape.Rmd](https://github.com/mcaitlinv/cawa-breeding/blob/main/05_au/scripts/genoscape.Rmd)
  
After identifying K=3 as the best K value for the adaptive loci, I used genoscapeRTools and tess3 to map the Q values from ADMIXTURE onto the map of the breeding range. This color coded K groups and used transparency to indicate where probabilities of assignment were less likely to produce a map of the likely adaptive groups across the breeding range. From these groups I created shapefiles of each group's outline to use as rasters for future analyses. 
