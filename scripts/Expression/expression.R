library(applyBy)
library(reshape2)

pop <- read.csv("capsule/code/CCanalysis/scripts//Expression/lineage_specifity.csv")
#pop <-  subset(pop, Lineage %in% "B")
#pop$ass_pop <-"Bcells"
sub_imm <- read.csv("capsule/code/CCanalysis/scripts/Expression/ImmExpBM.csv",header=T)
sub_imm <- t(sub_imm)
colnames(sub_imm) <- sub_imm[1,]
sub_imm <- sub_imm[-1,]
sub_imm <- apply(sub_imm,c(1,2),function(x){x<-as.numeric(x)})
rownames(sub_imm) <- sapply(rownames(sub_imm),function(x){gsub("X","",x)})

sub_imm <- apply(sub_imm,c(1,2),function(x){
  if(x>log2(47)){
    x<-1
  }
  else
  {
    x<-0
  }
})

sub_imm <- t(sub_imm)
sub_imm <- as.data.frame(sub_imm)
sub_imm$"pop"<- sapply(rownames(sub_imm),function(x){pop$ass_pop[match(x,pop$Cell)]})
sub_imm <- colSumsBy(as.matrix(sub_imm[,1:21755]), sub_imm$pop)
sub_imm <- t(sub_imm)
#sub_imm <- as.data.frame(sub_imm[which(rowSums(sub_imm)>0),])

sub_imm <- as.data.frame(sub_imm)
sub_imm$gene <- rownames(sub_imm)

#sub_imm$phen <- "Bcells"

sub_imm <- sub_imm %>%
              select(gene, phen, Bcells)

write.csv(sub_imm, "capsule/code/CCanalysis/scripts/Expression/expression_long_format_Bcells.csv")


