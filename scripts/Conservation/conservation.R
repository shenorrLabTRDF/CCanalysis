rm(list = ls())
require(ggplot2)
options(stringsAsFactors = FALSE)

genes <- read.csv("CCanalysis/scripts/Conservation//knownGeneOld9.csv")
genes <- genes[,1:5]
genes$X1 <- sapply(genes$X1,function(x){strsplit(x,".",fixed=TRUE)[[1]][1]})
genestoEnsm <- read.table("CCanalysis/scripts/Conservation//knownToEnsembl.txt")
genestoEnsm$V1 <- sapply(genestoEnsm$V1,function(x){strsplit(x,".",fixed=TRUE)[[1]][1]})

myGenes <- read.table("CCanalysis/scripts/Conservation/mart_export_all_genes.txt",sep=",")
colnames(myGenes) <- myGenes[1,]
myGenes<-myGenes[-1,]
myGenes<-as.data.frame(myGenes)
myGenes$"weirdName"<-sapply(myGenes$`Transcript stable ID`,function(x){genestoEnsm$V1[match(x,genestoEnsm$V2)]})
myGenes<-myGenes[complete.cases(myGenes),]

Cons60<-read.table("CCanalysis/scripts/Conservation/phastCons60way.txt",sep="\t")
myGenes<-cbind(myGenes,ncol=1)
colnames(myGenes)[4]<-"avScore"
myGenes<-myGenes[which(myGenes$weirdName %in% genes$X1),]

for (r in 1:nrow(myGenes)){
  g<-myGenes$weirdName[r]
  up<-genes[which(genes$X1 %in% g),4]
  to<-genes[which(genes$X1 %in% g),5]
  chr<-genes[which(genes$X1 %in% g),2]
  subCons1<-Cons60[which(Cons60[,2] %in% chr),]
  subCons1<-subCons1[which(subCons1$V3>=up),]
  subCons1<-subCons1[which(subCons1$V3<=to),] 
  
  subCons2<-Cons60[which(Cons60[,2] %in% chr),]
  subCons2<-subCons2[which(subCons2$V4>=up),]
  subCons2<-subCons2[which(subCons2$V4<=to),]
  
  subCons3<-Cons60[which(Cons60[,2] %in% chr),]
  subCons3<-subCons3[which(subCons3$V3<=up),]
  subCons3<-subCons3[which(subCons3$V4>=to),]

  
  subCons=rbind(subCons1,subCons2,subCons3)
  
  
  myGenes$avScore[r]<-mean(subCons$V13)
}

write.csv(myGenes,"CCanalysis/scripts/Conservation/mouse_conservation_scores_all.csv")

