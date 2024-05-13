rm(list = ls())
require(ggplot2)
options(stringsAsFactors = FALSE)


genes<-read.csv("capsule/code/CCanalysis/scripts/HPA_data/knownGeneOld11.csv",header=F)
genes<-genes[,1:5]
colnames(genes)<-genes[1,]
genes<-genes[-1,]
colnames(genes)<-c("Transcript stable ID version","chr","directionality","start","end")
myGenesVal<-read.table("capsule/code/CCanalysis/scripts/HPA_data/mart_export_human_genes.txt",sep=",")
colnames(myGenesVal)<-myGenesVal[1,]
myGenesVal<-myGenesVal[-1,]
myGenesVal<-as.data.frame(myGenesVal)
Cons60<-read.delim2("capsule/code/CCanalysis/scripts/HPA_data/phastCons100way.txt", header = FALSE)
myGenesVal<-left_join(myGenesVal,genes,by=c("Transcript stable ID version"))
myGenesVal<-myGenesVal[which(myGenesVal$`Gene name` %in% orru_study$Gene.annotation....1Mb.),]
myGenesVal<-myGenesVal[complete.cases(myGenesVal),]
myGenesVal<-cbind(myGenesVal,ncol=1)
colnames(myGenesVal)[10]<-"avScore"

for (r in 1:nrow(myGenesVal))
{
  up <- myGenesVal$start[r]
  to <- myGenesVal$end[r]
  subCons1 <- Cons60[which(Cons60[,2] %in% myGenesVal$chr[r]),]
  subCons1 <- subCons1[which(subCons1$V3>=up),]
  subCons1<-subCons1[which(subCons1$V3<=to),] 
  
  subCons2<-Cons60[which(Cons60[,2] %in% myGenesVal$chr[r]),]
  subCons2<-subCons2[which(subCons2$V4>=up),]
  subCons2<-subCons2[which(subCons2$V4<=to),]
  
  subCons3<-Cons60[which(Cons60[,2] %in% myGenesVal$chr[r]),]
  subCons3<-subCons3[which(subCons3$V3<=up),]
  subCons3<-subCons3[which(subCons3$V4>=to),]
  
  
  subCons=rbind(subCons1,subCons2,subCons3)
  myGenesVal$avScore[r]<-mean(as.numeric(subCons$V13))
}

write.csv(myGenesVal,"CC/scripts/HPA_data/orru_study_conservation.csv")


