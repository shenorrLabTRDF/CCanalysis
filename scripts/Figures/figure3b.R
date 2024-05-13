library("reshape2")
library("ggplot2")
install.packages("doBy")
require("doBy")
library("dplyr")
library("plyr")

#mboat1

conservation_analysis <- read.csv("capsule/code/CCanalysis/scripts/SignalPropagation/conservation_analysis_updated_adj_005.csv")
trans <- unique(conservation_analysis$Gene.name[which(conservation_analysis$tr_expr %in% "trans")])
trans <- c(trans,"2610019F03Rik")
trans <- "Arhgef37"
pop <- read.csv("capsule/code/CCanalysis/scripts/Expression/lineage_specifity.csv")
sub_imm <- read.csv("capsule/code/CCanalysis/scripts/Expression/ImmExpBM.csv",header=T)
sub_imm <- t(sub_imm)
colnames(sub_imm) <- sub_imm[1,]
sub_imm <- sub_imm[-1,]
sub_imm <- apply(sub_imm,c(1,2),function(x){x<-as.numeric(x)})
rownames(sub_imm) <- sapply(rownames(sub_imm),function(x){gsub("X","",x)})

sub_imm <- subset(sub_imm, rownames(sub_imm) %in% trans)

sub_imm <- t(sub_imm)
sub_imm <- as.data.frame(sub_imm)
sub_imm$"pop"<- sapply(rownames(sub_imm),function(x){pop$ass_pop[match(x,pop$Cell)]})
sub_imm<- melt(sub_imm)

expr_max <- summaryBy(value~pop+variable, data = sub_imm, FUN = max)
expr_max$"expr"<-0
expr_max$expr[which(expr_max$value.max>log2(47))]<-1

pop <- read.csv("capsule/code/CCanalysis/scripts/Expression/lineage_specifity.csv")
pop <-  subset(pop, Lineage %in% "B")
pop$ass_pop <-"Bcells"

sub_imm <- read.csv("capsule/code/CCanalysis/scripts/Expression/ImmExpBM.csv",header=T)
sub_imm <- t(sub_imm)
colnames(sub_imm) <- sub_imm[1,]
sub_imm <- sub_imm[-1,]
sub_imm <- apply(sub_imm,c(1,2),function(x){x<-as.numeric(x)})
rownames(sub_imm) <- sapply(rownames(sub_imm),function(x){gsub("X","",x)})
sub_imm <- subset(sub_imm, rownames(sub_imm) %in% trans)

sub_imm <- t(sub_imm)
sub_imm <- as.data.frame(sub_imm)
sub_imm$"pop"<- sapply(rownames(sub_imm),function(x){pop$ass_pop[match(x,pop$Cell)]})
sub_imm<- sub_imm[complete.cases(sub_imm),]
sub_imm<- melt(sub_imm)

expr_maxB <- summaryBy(value~pop+variable, data = sub_imm, FUN = max)
expr_maxB$"expr"<-0
expr_maxB$expr[which(expr_maxB$value.max>log2(47))]<-1

expr_max_bind <- rbind(expr_max,expr_maxB)
expr_max_bind$pop[which(expr_max_bind$pop %in% "restB")] <- "lateB"
expr_max_bind$pop[which(expr_max_bind$pop %in% "ly6gPos")] <- "Granulocytes"
expr_max_bind$pop[which(expr_max_bind$pop %in% "classicalMonocytes")] <- "Monocytes"
expr_max_bind$pop[which(expr_max_bind$pop %in% "lassicalMonocytes")] <- "Monocytes"
expr_max_bind <- expr_max_bind %>%
  filter(!pop %in% c("NKcells","preproB"))
expr_max_bind$pop <- factor(x= expr_max_bind$pop, levels =c("StemCells ","proB","Monocytes","lateB","Granulocytes","CD8Tcells","CD4Tcells","Bcells"))

expr_max_bind <- expr_max_bind %>%
  mutate(pop = fct_relevel(pop, "StemCells ","proB","Monocytes","lateB","Granulocytes","CD8Tcells","CD4Tcells","Bcells"))

ggplot(expr_max_bind,aes(x=pop, y=value.max))  + 
  geom_bar(position="dodge",stat="identity")+ 
  geom_hline(yintercept = log2(47)) +
  #scale_y_continuous(breaks = levels(seq(0,8,by=2))) +
  coord_flip()
  #facet_wrap(~variable)  

###3b right
conservation_analysis_updated_adj_005 <- read.csv("capsule/code/CCanalysis/scripts/SignalPropagation/conservation_analysis_updated_adj_005.csv")
cis<- unique(conservation_analysis_updated_adj_005$Gene.name[which(conservation_analysis_updated_adj_005$tr_expr %in% "cis")])
cis <- cis[100:119]
cis <- c(cis)

all_associations <- read.csv("capsule/code/CCanalysis/scripts/Figures/all2.csv")
all_associations <- all_associations %>%
  dplyr::filter(!phen %in% c("preproB","NKcells"))
###ading extra snps tp cover low associated snps
prob<- read.csv("capsule/code/CCanalysis/scripts/mapping/results/ProB.csv")
prob$"phen"<-"proB"
cd8t<- read.csv("capsule/code/CCanalysis/scripts/mapping/results/CD8Tcells.csv")
cd8t$"phen" <- "CD8Tcells"

all_associations <- rbind(all_associations,prob,cd8t)
all_associations <- all_associations[which(!duplicated(all_associations)),]

imm_exons_intervals_rel <- read.csv("capsule/code/CCanalysis/scripts/Loci_filtering/rel_intervals_updated.csv")
all_associations <- subset(all_associations, marker %in% imm_exons_intervals_rel$marker)

intervals <- subset(imm_exons_intervals_rel, Gene.name %in% conservation_analysis_updated_adj_005$Gene.name)
intervals$"adj" <- paste0(intervals$marker,"_",intervals$Gene.name)
intervals <- subset(intervals,!duplicated(intervals$adj))

all_associations <- subset(all_associations, marker %in% intervals$marker)
all_associations <- left_join(all_associations,intervals, by="marker",relationship="many-to-many")
all_associations_agg <- all_associations %>%
                            dplyr::group_by(marker) %>%
                            dplyr::summarise(n_phen = n_distinct(phen))
all_associations <- subset(all_associations, Gene.name %in% trans)
#if there is more than one snp, chose on with higher lod
all_associations <- ddply(all_associations,.(phen,Gene.name), function(markers){
  if(nrow(markers)>1){
    markers[which(markers$lod %in% max(markers$lod)),]
  }else{
    markers
  }
})

all_associations <- all_associations %>%
                     dplyr::select(phen,lod,Gene.name)
phens <- unique(all_associations$phen)

fdr_key <- data.frame(cell_types = c("Bcells","CD4Tcells","CD8Tcells","Granulocytes","lateB","Monocytes","NKcells","proB","StemCells"),
                      fdr_tr = c(5.75,7,8.5,6,6,6,8.9,6.6,10.11)) 

all_associations$"tr" <- sapply(all_associations$phen,function(x){fdr_key$fdr_tr[match(x,fdr_key$cell_types)]})
all_associations$phen <- factor(x= all_associations$phen, levels =c("StemCells","proB","Monocytes","lateB","Granulocytes","CD8Tcells","CD4Tcells","Bcells"))


# to_add <- data.frame(phens[which(!phens %in% all_associations$phen[which(all_associations$Gene.name %in% "Mboat1")])],0.01,"Mboat1")
# colnames(to_add) <- colnames(all_associations)
# all_associations <- rbind(all_associations,to_add)
ggplot(all_associations)+
  geom_bar(aes(x=phen, y=lod),stat="identity")+ 
  geom_line(aes(x=phen, y=tr,group=1),stat="identity",color="blue",size=0.5)+
  coord_flip()
  facet_wrap(~Gene.name)
