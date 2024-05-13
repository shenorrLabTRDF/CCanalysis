library("dplyr")
library("plyr")
library("ggplot2")
install.packages("readxl")
library("readxl")
options(stringsAsFactors = FALSE)

rel_genes <- read.csv("capsule/code/CCanalysis/scripts/Loci_filtering/rel_genes_remapped.csv")
rel_genes$"val_pvalue" <- 0
colnames(rel_genes)[14] <- 'val_pvalue'


mapfile <- read.csv("capsule/code/CCanalysis//scripts/validation/mapfile.csv")
hap <- readRDS("capsule/code/CCanalysis/scripts/Loci_filtering/hap_reconstruction_do.RData")
Supplementary_table_1_Mouse_Strains <- read_excel("capsule/code/CCanalysis/scripts/validation/Supplementary table 1- Mouse Strains.xlsx")


####prepare phen
phen <- read.csv("capsule/code/CCanalysis//scripts/validation/phen_2dataset_old_gating_without_overlappingStrains.csv")

###filter for immune system genes and cc variation
strict <- read.csv("capsule/code/CCanalysis/scripts/validation/strict_matrix_second_set_noOverlap.csv")
phen <- subset(phen,! alias %in% Supplementary_table_1_Mouse_Strains$Alias[which(Supplementary_table_1_Mouse_Strains$dataset1==1 & Supplementary_table_1_Mouse_Strains$dataset2==1)])
strict <- subset(strict,! Alias %in% Supplementary_table_1_Mouse_Strains$Alias[which(Supplementary_table_1_Mouse_Strains$dataset1==1 & Supplementary_table_1_Mouse_Strains$dataset2==1)])


for (r in 1:nrow(rel_genes)){
  sub_phen <- cbind(as.character(phen$alias),phen[,as.character(rel_genes$phen[r])])
  sub_phen <- as.data.frame(sub_phen)
  colnames(sub_phen) <- c("strain","phen")
  ##1.subet phen using strict matrix for specific phenotype
  sub_phen <- subset(sub_phen,sub_phen$strain %in% strict$Alias[strict[,as.character(rel_genes$phen[r])]!=0])
  ###add founder column
    sub_phen$'founder'<- sapply(sub_phen$strain,function(x){hap[match(x,hap$X),rel_genes$marker[r]]})
    ##added mean calculation step
    sub_phen$phen <- as.numeric(sub_phen$phen)
    ###add p.value to initial matrix
    results <- summary(lm(sub_phen$phen~as.factor(sub_phen$founder)))
    f <- results$fstatistic
    p_value <- pf(f[1],f[2],f[3],lower.tail=FALSE)
    rel_genes$'val_pvalue'[r] <- p_value
}

rel_genes$"adjusted" <- p.adjust(rel_genes$val_pvalue,method =  "BH")
rel_genes2 <- subset(rel_genes, adjusted < 0.05)

# rel_genes <- rel_genes %>%
#   dplyr::mutate(adj_pvalue =  p.adjust(val_pvalue,method =  "BH")) %>%
#   dplyr::filter(val_pvalue < 0.05)
write.csv(rel_genes,"capsule/code/CCanalysis/scripts/validation/rel_genes_validated_updated.csv")
write.csv(rel_genes2,"capsule/code/CCanalysis/scripts/validation/rel_genes_validated_updated_validated_005.csv")

##creating mean matrix for first cohort as well

phen1 <- read.csv("capsule/code/CCanalysis/scripts/Loci_filtering/phen_BM_updated_20200711.csv")
ms_pheno <- c(11,38)
phen1 <- phen1[-ms_pheno,]
colnames(phen1)[1]<-"alias"
colnames(phen1)[7]<-"proB"
colnames(phen1)[8]<-"earlyB"
colnames(phen1)[9]<-"lateB"

strict1 <- read.csv("capsule/code/CCanalysis/scripts/Loci_filtering/strict_matrix_old.csv")
strict1 <- strict1 %>%
  dplyr::select(samples,Tcells,NKcells,ly6gPos,classicalMonocytes,proB,earlyB,lateB,CD8Tcells,CD4Tcells,StemCells,Bcells)

colnames(strict1)[2:12] <- colnames(phen1)[3:13]
colnames(strict1)[1] <- "Alias"

 phen1 <- subset(phen1,! alias %in% Supplementary_table_1_Mouse_Strains$Alias[which(Supplementary_table_1_Mouse_Strains$dataset1==1 & Supplementary_table_1_Mouse_Strains$dataset2==1)])
 strict1 <- subset(strict1,! Alias %in% Supplementary_table_1_Mouse_Strains$Alias[which(Supplementary_table_1_Mouse_Strains$dataset1==1 & Supplementary_table_1_Mouse_Strains$dataset2==1)])
 strict1 <- subset(strict1, !duplicated(Alias)) 

founder_avg1 <- matrix(ncol=8,nrow=0)
for (r in 1:nrow(rel_genes2))
{
  sub_phen <- cbind(as.character(phen1$alias),phen1[,as.character(rel_genes2$phen[r])])
  sub_phen <- as.data.frame(sub_phen)
  colnames(sub_phen) <- c("strain","phen")
  ##1.subet phen using strict matrix for specific phenotype
  sub_phen <- subset(sub_phen,sub_phen$strain %in% strict1$Alias[strict1[,as.character(rel_genes2$phen[r])]!=0])
  ###add founder column
  sub_phen$'founder'<- sapply(sub_phen$strain,function(x){hap[match(x,rownames(hap)),rel_genes2$marker[r]]})
  ##added mean calculation step
  sub_phen$phen <- as.numeric(sub_phen$phen)
  temp <- sub_phen %>%
    dplyr::group_by(founder) %>%
    dplyr::summarize(avg=mean(phen), na.rm = TRUE) %>%
    dplyr::mutate(snp=rel_genes2$marker[r]) %>%
    dplyr::mutate(phenotype=rel_genes2$phen[r])
  founder_avg1 <- rbind(founder_avg1, temp)
}

colnames(founder_avg1)[2] <- "avg_first"

founder_avg <- matrix(ncol=8,nrow=0)
for (r in 1:nrow(rel_genes2)){
  sub_phen <- cbind(as.character(phen$alias),phen[,as.character(rel_genes2$phen[r])])
  sub_phen <- as.data.frame(sub_phen)
  colnames(sub_phen) <- c("strain","phen")
  ##1.subet phen using strict matrix for specific phenotype
  sub_phen <- subset(sub_phen,sub_phen$strain %in% strict$Alias[strict[,as.character(rel_genes2$phen[r])]!=0])
  ###add founder column
  sub_phen$'founder'<- sapply(sub_phen$strain,function(x){hap[match(x,rownames(hap)),rel_genes2$marker[r]]})
  ##added mean calculation step
  sub_phen$phen <- as.numeric(sub_phen$phen)
  temp <- sub_phen %>%
    dplyr::group_by(founder) %>%
    dplyr::summarize(avg=mean(phen), na.rm = TRUE) %>%
    dplyr::mutate(snp=rel_genes2$marker[r]) %>%
    dplyr::mutate(phenotype=rel_genes2$phen[r])
  founder_avg <- rbind(founder_avg, temp)
}

founder_avg <- founder_avg %>%
                 dplyr::select(-na.rm)
founder_avg1 <- founder_avg1 %>%
  dplyr::select(-na.rm)

founders <- full_join(founder_avg, founder_avg1,relationship = "many-to-many")
founders <- founders[complete.cases(founders),]


founders_cor <- ddply(founders, .(snp,phenotype), function(ass){
  if(nrow(ass) > 4){
  cor(ass$avg_first,ass$avg,method = "spearman")}
})

colnames(founders_cor)[3] <- "correlation"
#founders_cor <- subset(founders_cor, !phenotype %in% "StemCells")
ggplot(founders_cor, aes(x=correlation,fill=phenotype))+geom_density()


sub_0 <- subset(founders_cor, correlation > 0)
sub_0.5 <- subset(founders_cor, correlation > 0.5)

genes0 <- subset(rel_genes2, phen %in% sub_0$phenotype & marker %in% sub_0$snp)
genes05 <- subset(rel_genes, phen %in% sub_0.5$phenotype & SNP %in% sub_0.5$snp)
write.csv(genes0,"capsule/code/CCanalysis/scripts/validation/gene_cor_greater_0_updated.csv")
write.csv(genes05,"CC/scripts/validation/gene_cor_greater_05.csv")
