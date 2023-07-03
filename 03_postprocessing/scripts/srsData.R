library(remotes)
remotes::install_github("eriqande/srsStuff", upgrade = FALSE)
library(srsStuff)
library(tidyverse)


snames <- read_lines("cawa.all.merged.gatk.filtered.miss0.75.recode.vcf.nohead")
gc()
allCAWA<-srs_covar(file="allele_depth_0.75.txt", sample_names=snames, freq_thresh = 0)
saveRDS(allCAWA, "cawa.srs.filtered.rds")

