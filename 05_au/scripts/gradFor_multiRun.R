library(extendedForest)
library(gradientForest)
require(gstat)
require(RColorBrewer)
library(tidyverse)
library(readr)
require(psych)

##Make 5 different sets of 50k random SNPs
myargs <- commandArgs(trailingOnly=T)
x <- myargs[1]
for (val in x) {
  name<-paste("cawa.16ClimGroup.minInd0.7.majmin4.scaf.noNA.50Krandom", val,".sort", sep="")
  # gen<-read_delim(paste("results/au/gradientFor/genomeData/", name,".txt", sep=""),delim="\t",col_names = T)
  # Gcawa<-t(gen)
  # write.table(Gcawa,paste("results/au/gradientFor/genomeData/", name, ".wide.csv", sep=""),row.names=F,quote=F,sep=",")
  Gcawa<- read_delim(paste("results/au/gradientFor/genomeData/", name, ".wide.csv", sep=""), delim=",")
  ##Using each of the sets of 50k random SNPs, make gradientforest model
  Ecawa <- read_delim("results/au/gradientFor/cawaEnv.220627.txt",delim="\t", col_names=T)
  ##RE-do name to include date
  name<-paste("cawa.16ClimGroup.minInd0.7.majmin4.scaf.noNA.50Krandom", val,".sort.220627", sep="")
  ##Remove lat/long and qscat
  Ecawa <- Ecawa %>% select(!c( lat, long, qscat))
  # ##Specify predictors which are the environmentals vars
  # preds <- colnames(Ecawa)
  # print(preds)
  # specs <- colnames(Gcawa)
  
  # nSites <- as.numeric(dim(Gcawa)[1])
  # nSpecs <- as.numeric(dim(Gcawa)[2])
  
  # # set depth of conditional permutation
  # lev <- floor(log2(nSites*0.368/2))
  
  
  # ##Making the gradient forest model with env as pred and response as genomic sites
  # ##This takes awhile and spits out thousands of lines of errors
  # ##Ignore the errors
  # cawaforest=gradientForest(cbind(Ecawa,Gcawa), predictor.vars=preds, response.vars=specs, ntree=500, transform = NULL, compact=T,nbin=101, maxLevel=lev,trace=T)
  
  # saveRDS(cawaforest, file = paste("results/au/gradientFor/cawaforest", name,".rds", sep=""))
  cawaforest<-readRDS(paste("results/au/gradientFor/cawaforest", name,".rds", sep=""))
  #Check predictor importance
  # allpredictors=names(importance(cawaforest))
  # ##Save as pdf 
  # pdf(paste("results/au/gradientFor/output/", name, ".importance.pdf", sep=""))
  # plot(cawaforest,plot.type="O")
  # dev.off()
  
  # ##Check for correlation in variables
  # ##Aiming to take the top 10 variables and pare down until max number of variables with corr <0.7
  # cor<-cor(Ecawa,method="pearson")
  # cor %>% write.table(paste("results/au/gradientFor/tables/", name, ".corr", sep=""),quote=F,sep="\t",row.names=T)
  
  ##Start randomization to check how accurate model is
  realtotal=cawaforest$species.pos.rsq
  realaverage=sum(cawaforest$result)/realtotal
  
  
  #Create randomization values for 100 different runs of random values
  #Paste into table
  for (i in 1:100) {
    ##Sample randomly the rows of the environment
    predsR=Ecawa[sample(nrow(Ecawa)),]
  
    ##Build your model with the randomly sampled rows of the environment
    ##These will be mismatched to the rows of genomic data, thereby randomizing your predictors
    cawaforestR=gradientForest(cbind(predsR,Gcawa), predictor.vars=colnames(predsR), response.vars=colnames(Gcawa), ntree=100, transform = NULL, compact=F, maxLevel=1,trace=F)
    ##Find the random total r-sq and average r-sq
    randtotal=cawaforestR$species.pos.rsq
    randaverage=sum(cawaforestR$result)/randtotal
  
    ##Write random values into tables
    write.table(randtotal,file=paste("results/au/gradientFor/tables/",name,"bslmm.nobin.randtotal",sep=""),row.names=FALSE,col.names=FALSE,append=T)
    write.table(randaverage,file=paste("results/au/gradientFor/tables/",name,"bslmm.nobin.randaverage",sep=""),row.names=FALSE,col.names=FALSE,append=T)
  }
  
  ##Read in the 500 random averages and r-squareds
  randa=read.table(paste("results/au/gradientFor/tables/",name,"bslmm.nobin.randtotal",sep=""))
  randt=read.table(paste("results/au/gradientFor/tables/",name,"bslmm.nobin.randaverage",sep=""))
 
  quantile(randa$V1,c(0.995))
  
  hist(randa$V1,main="Average R-squared of Random Gradient Forests",xlab="Average R-squared of SNPs", xlim =c(0.09, 0.18))

  abline(v=realaverage,col="red")
  ##Make a graph showing the distribution of randoms
  pdf(paste("results/au/gradientFor/output/", name, ".randomization.hist.par.pdf", sep=""))
  par(mfrow=c(2,1))
  hist(randa$V1,main="Average R-squared of Random Gradient Forests",xlab="Average R-squared of SNPs", xlim =c(0.09, 0.18))
  abline(v=realaverage,col="red") #replace with your actual number from realaverage
  dev.off()
 }
