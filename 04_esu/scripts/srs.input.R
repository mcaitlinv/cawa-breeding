#remotes::install_github("eriqande/srsStuff", upgrade = FALSE)
library(srsStuff)
library(tidyverse)


snames <- read_lines("samples1X.list")
CAWA.filt<-srs_covar(file="copy/cawa.bqsr.miss0.2.qual30.1X.fst85.filtRda.recode.AD.txt", sample_names=snames, freq_thresh = 0)
saveRDS(CAWA.filt, "copy/cawa.srs.filt.85fst.RDA.rds")
