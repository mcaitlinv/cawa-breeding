rm(list=ls())
require(extendedForest)
require(gradientForest)
library(tidyverse)
require(readr)


gen<-read_delim("results/au/gradientForAdaptOnly/cawa.adaptOnly.16ClimGroup.minInd0.5.majmin4.scaf.noNA.txt",delim="\t",col_names = F)
gen[is.na(gen)] <- 0
gen <- gen[, -1:-2]
Gcawa<-t(gen)
write.table(Gcawa,"results/au/gradientForAdaptOnly/genomeData/cawa.adaptOnly.16ClimGroup.GcawaminInd0.5.majmin4.noscaf.noNA.wide.csv",row.names=F,quote=F,sep=",")
Gcawa<- read_delim("results/au/gradientForAdaptOnly/genomeData/cawa.adaptOnly.16ClimGroup.GcawaminInd0.5.majmin4.noscaf.noNA.wide.csv", delim=",")
Ecawa <- read_delim("results/au/gradientFor/cawaEnv.220627.txt",delim="\t", col_names=T)
##Remove lat/long and qscat and 
Ecawa <- Ecawa %>% select(bio10,bio13,bio15,tree)
##Specify predictors which are the environmentals vars
preds <- colnames(Ecawa)
print(preds)
specs <- colnames(Gcawa)
  
nSites <- as.numeric(dim(Gcawa)[1])
nSpecs <- as.numeric(dim(Gcawa)[2])
  
# set depth of conditional permutation
lev <- floor(log2(nSites*0.368/2))
 
##Making the gradient forest model with env as pred and response as genomic sites
##This takes awhile and spits out thousands of lines of errors
##Ignore the errors
#cawaforest=gradientForest(cbind(Ecawa,Gcawa), predictor.vars=preds, response.vars=specs, ntree=500, transform = NULL, compact=T,nbin=101, maxLevel=lev,trace=T)
#saveRDS(cawaforest, file = "results/au/gradientForAdaptOnly/cawaAdaptForest.rds")

cawaforest<-readRDS("results/au/gradientForAdaptOnly/cawaAdaptForest.rds")
#Check predictor importance
allpredictors=names(importance(cawaforest))
##Save as pdf 
pdf("results/au/gradientForAdaptOnly/output/cawaAdaptForest.importance.pdf")
plot(cawaforest,plot.type="O")
dev.off()
  
##Check for correlation in variables
##Aiming to take the top 10 variables and pare down until max number of variables with corr <0.7
cor<-cor(Ecawa,method="pearson")
write.table(cor, "results/au/gradientForAdaptOnly/tables/cawaAdaptForest.corr",quote=F,sep="\t",row.names=T)
  
##Start randomization to check how accurate model is
realtotal=cawaforest$species.pos.rsq
realaverage=sum(cawaforest$result)/realtotal
realaverage

##Create randomization values for 100 different runs of random values
##Paste into table
for (i in 1:100){

##Sample randomly the rows of the environment
predsR=Ecawa[sample(nrow(Ecawa)),]

##Build your model with the randomly sampled rows of the environment
##These will be mismatched to the rows of genomic data, thereby randomizing your predictors
cawaforestR=gradientForest(cbind(predsR,Gcawa), predictor.vars=colnames(predsR), response.vars=colnames(Gcawa), ntree=100, transform = NULL, compact=F, maxLevel=1,trace=F)

##Find the random total r-sq and average r-sq
randtotal=cawaforestR$species.pos.rsq
randaverage=sum(cawaforestR$result)/randtotal

##Write random values into tables
write.table(randtotal,file="results/au/gradientForAdaptOnly/tables/cawa.adaptOnly.bslmm.nobin.randtotal",row.names=FALSE,col.names=FALSE,append=T)
write.table(randaverage,file="results/au/gradientForAdaptOnly/tables/cawa.adaptOnly.bslmm.nobin.randaverage",row.names=FALSE,col.names=FALSE,append=T)
}

##Read in the 100 random averages and r-squareds
randt=read.table("results/au/gradientForAdaptOnly/tables/cawa.adaptOnly.bslmm.nobin.randtotal",sep="")
randa=read.table("results/au/gradientForAdaptOnly/tables/cawa.adaptOnly.bslmm.nobin.randaverage",sep="")
quantile(randa$V1,c(0.995))

##Make a graph showing the distribution of randoms
pdf("results/au/gradientForAdaptOnly/output/cawaAdaptForest.randomization.100trees.nobin.bslmm.hist.par.pdf")
par(mfrow=c(2,1))
hist(randa$V1,main="Average R-squared of Random Gradient Forests",xlab="Average R-squared of SNPs", xlim =c(0.09, 0.18))
abline(v=realaverage,col="red") #replace with your actual number from realaverage
dev.off()

realaverage> quantile(randa$V1,c(0.995))
