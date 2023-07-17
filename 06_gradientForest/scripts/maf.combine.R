rm(list=ls())
require(readr)
library(tidyverse)

ABFMK<-read_delim("results/au/gradientForAdaptOnly/AB.FortMacKay.minInd0.5.majmin4.adapt.mafs.gz",delim="\t") %>% rename(ABFMK_freq=knownEM) %>% select(chromo,position,major,ref,ABFMK_freq)

ABFMM<-read_delim("results/au/gradientForAdaptOnly/AB.FortMcMurray.minInd0.5.majmin4.adapt.mafs.gz",delim="\t" )%>% rename(ABFMM_freq=knownEM) %>% select(chromo,position,major,ref,ABFMM_freq)

ABSL <-read_delim("results/au/gradientForAdaptOnly/AB.SlaveLake.minInd0.5.majmin4.adapt.mafs.gz",delim="\t" )%>% rename(ABSL_freq=knownEM) %>% select(chromo,position,major,ref,ABSL_freq)

MBR <-read_delim("results/au/gradientForAdaptOnly/MB.Rennie.minInd0.5.majmin4.adapt.mafs.gz",delim="\t" )%>% rename(MBR_freq=knownEM) %>% select(chromo,position,major,ref,MBR_freq)

MNF <-read_delim("results/au/gradientForAdaptOnly/MN.Finland.minInd0.5.majmin4.adapt.mafs.gz",delim="\t" )%>% rename(MNF_freq=knownEM) %>% select(chromo,position,major,ref,MNF_freq)

NBMA <-read_delim("results/au/gradientForAdaptOnly/NB.McAdam.minInd0.5.majmin4.adapt.mafs.gz",delim="\t") %>%  rename(NBMA_freq=knownEM) %>% select(chromo,position,major,ref,NBMA_freq)

NCO <-read_delim("results/au/gradientForAdaptOnly/NC.Otto.minInd0.5.majmin4.adapt.mafs.gz",delim="\t" )%>% rename(NCO_freq=knownEM) %>% select(chromo,position,major,ref,NCO_freq)

NHC <-read_delim("results/au/gradientForAdaptOnly/NH.Canaan.minInd0.5.majmin4.adapt.mafs.gz",delim="\t" )%>% rename(NHC_freq=knownEM) %>% select(chromo,position,major,ref,NHC_freq)

NYMR <-read_delim("results/au/gradientForAdaptOnly/NY.MooseRiver.minInd0.5.majmin4.adapt.mafs.gz",delim="\t" )%>% rename(NYMR_freq=knownEM) %>% select(chromo,position,major,ref,NYMR_freq)

PABR <-read_delim("results/au/gradientForAdaptOnly/PA.BrightRun.minInd0.5.majmin4.adapt.mafs.gz",delim="\t" )%>% rename(PABR_freq=knownEM) %>% select(chromo,position,major,ref,PABR_freq)

PAWR <-read_delim("results/au/gradientForAdaptOnly/PA.WolfRun.minInd0.5.majmin4.adapt.mafs.gz",delim="\t" )%>% rename(PAWR_freq=knownEM) %>% select(chromo,position,major,ref,PAWR_freq)

QCLT <-read_delim("results/au/gradientForAdaptOnly/QC.Laterriere.minInd0.5.majmin4.adapt.mafs.gz",delim="\t" )%>% rename(QCLT_freq=knownEM) %>% select(chromo,position,major,ref,QCLT_freq)

QCSF <-read_delim("results/au/gradientForAdaptOnly/QC.St-Fulgence.minInd0.5.majmin4.adapt.mafs.gz",delim="\t" )%>% rename(QCSF_freq=knownEM) %>% select(chromo,position,major,ref,QCSF_freq)

RISF <-read_delim("results/au/gradientForAdaptOnly/RI.SpragueFarm.minInd0.5.majmin4.adapt.mafs.gz",delim="\t" )%>% rename(RISF_freq=knownEM) %>% select(chromo,position,major,ref,RISF_freq)

WIEC <-read_delim("results/au/gradientForAdaptOnly/WI.EauClaire.minInd0.5.majmin4.adapt.mafs.gz",delim="\t" )%>% rename(WIEC_freq=knownEM) %>% select(chromo,position,major,ref,WIEC_freq)

WIRW <-read_delim("results/au/gradientForAdaptOnly/WV.Richwood.minInd0.5.majmin4.adapt.mafs.gz",delim="\t" )%>% rename(WIRW_freq=knownEM) %>% select(chromo,position,major,ref,WIRW_freq)

maf<- ABFMK %>% left_join(ABFMM) %>% left_join(ABSL) %>% left_join(MBR) %>% left_join(MNF) %>% left_join(NBMA) %>% left_join(NCO) %>% left_join(NHC) %>% left_join(NYMR) %>% left_join(PABR) %>% left_join(PAWR) %>% left_join(QCLT) %>% left_join(QCSF) %>% left_join(RISF) %>% left_join(WIEC) %>% left_join(WIRW)

maf %>% select(-major,-ref) %>% write.table('results/au/gradientForAdaptOnly/cawa.adaptOnly.16ClimGroup.minInd0.5.majmin4.NA.txt',row.names=F,quote=F,sep="\t")

#maf <-read.table('/scratch/summit/caitlinv\@colostate.edu/cawa-breed-wglc/results/au/gradientFor/cawa.16ClimGroup.minInd0.5.majmin4.NA.txt', header=T)

maf[is.na(maf)] <- 0
maf2<- maf
maf2 %>% write.table('results/au/gradientForAdaptOnly/cawa.adaptOnly.16ClimGroup.minInd0.5.majmin4.freq.txt',row.names=F,quote=F,sep="\t")
# maf2<- read.table('results/au/gradientFor/cawa.16ClimGroup.minInd0.5.majmin4.freq.txt')

# maf2<- maf %>% select(-chromo,-position)
# maf3<- maf2[,apply(maf2, 2, var, na.rm=TRUE) !=0]
# dim(maf2)
# dim(maf3)
# dat_maf<-t(maf2)
# write.table(dat_maf,'results/au/gradientFor/cawa.16ClimGroup.minInd0.5.majmin4.allele.freq.csv',row.names=F,quote=F,sep=",")

#cat /scratch/summit/caitlinv\@colostate.edu/cawa-breed-wglc/results/au/gradientFor/cawa.16ClimGroup.minInd0.5.majmin4.NA.txt | cut -f 1-11 | sed 1d | awk 'BEGIN{FS="\t";OFS="\t"};{count=0;for(i=3; i<12; i++) {if($i!="NA") {count++}};if (count>=9){print $0}}' > /scratch/summit/caitlinv\@colostate.edu/cawa-breed-wglc/results/au/gradientFor/cawa.16ClimGroup.minInd0.5.majmin4.scaf.noNA.txt
