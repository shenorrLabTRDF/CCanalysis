library("dplyr")
library("plyr")
library("ggplot2")
install.packages("readxl")
library("readxl")
options(stringsAsFactors = FALSE)
library("cytoreason.cc.client")

rel_genes <- read.csv("capsule/code/CCanalysis/scripts/Loci_filtering/rel_genes_remapped.csv")
rel_genes$"val_pvalue" <- 0
colnames(rel_genes)[14] <- 'val_pvalue'


Supplementary_table_1_Mouse_Strains <- read_excel("capsule/code/CCanalysis/scripts/validation/Supplementary table 1- Mouse Strains.xlsx")

hap <- readRDS("capsule/code/CCanalysis/scripts/Loci_filtering/hap_reconstruction_do.RData")
rownames(hap) <-hap$X
hap <- hap[,-1,drop=FALSE]
hap <- hap[, colnames(hap) %in% rel_genes$marker]


####prepare phen
phen <- read.csv("capsule/code/CCanalysis//scripts/validation/phen_2dataset_old_gating_without_overlappingStrains.csv")

###filter for immune system genes and cc variation
strict <- read.csv("capsule/code/CCanalysis/scripts/validation/strict_matrix_second_set_noOverlap.csv")
phen <- subset(phen,! alias %in% Supplementary_table_1_Mouse_Strains$Alias[which(Supplementary_table_1_Mouse_Strains$dataset1==1 & Supplementary_table_1_Mouse_Strains$dataset2==1)])
strict <- subset(strict,! Alias %in% Supplementary_table_1_Mouse_Strains$Alias[which(Supplementary_table_1_Mouse_Strains$dataset1==1 & Supplementary_table_1_Mouse_Strains$dataset2==1)])

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

n_perm <- 100
res <- cytoreason.cc.client::sapply_dist(setNames(1:n_perm, 1:n_perm),
                                         FUN = function(n, phen, strict, phen1, strict1, rel_genes, hap){
  library("dplyr")
  library("plyr")
  
  row_names <- rownames(hap)
  #hap <- hap[sample(nrow(hap)), ]
  hap <- hap[sample(c(1:8),nrow(hap),replace=TRUE), ]
  rownames(hap) <- row_names
  # hap <- as.data.frame(t(hap))                                         
  # hap$gene <- sapply(rownames(hap),function(x){rel_genes$Gene.name[match(x,rel_genes$marker)]})    
  # hap$marker<- rownames(hap)
  # hap <- ddply(hap,.(gene), function(gen){
  #    gen_marker <- gen$marker
  #    mat <- do.call("rbind", replicate(nrow(gen), gen[1,,drop=FALSE], simplify = FALSE)) 
  #    mat$marker <- gen_marker
  #    mat
  # })
  # rownames(hap) <- hap$"marker"
  # hap <- hap %>% select(-c("gene","marker"))
  # hap <- as.data.frame(t(hap))  
  
  
  for (r in 1:nrow(rel_genes)){
    sub_phen <- cbind(as.character(phen$alias),phen[,as.character(rel_genes$phen[r])])
    sub_phen <- as.data.frame(sub_phen)
    colnames(sub_phen) <- c("strain","phen")
    ##1.subet phen using strict matrix for specific phenotype
    sub_phen <- subset(sub_phen,sub_phen$strain %in% strict$Alias[strict[,as.character(rel_genes$phen[r])]!=0])
    ###add founder column
    sub_phen$'founder'<- sapply(sub_phen$strain,function(x){hap[match(x,rownames(hap)),rel_genes$marker[r]]})
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
  if(nrow(rel_genes2)==0){
    return(rel_genes2)
  }else{  
  
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
 
  print(nrow(founder_avg))
  print(nrow(founder_avg1))
  
  if(nrow(founder_avg)>0 & nrow(founder_avg1)>0){

  founders <- full_join(founder_avg, founder_avg1,relationship = "many-to-many")
  founders <- founders[complete.cases(founders),]

  founders_cor <- ddply(founders, .(snp,phenotype), function(ass){
    if(nrow(ass) > 4){
      cor(ass$avg_first,ass$avg,method = "spearman")}
  })

  colnames(founders_cor)[3] <- "correlation"
  sub_0 <- subset(founders_cor, correlation > 0)
  # colnames(sub_0)[1:2] <- c("marker","phen")
  # rel_genes2 <- left_join(rel_genes2, sub_0)
  
  }else{
    sub_0<-data.frame(nrow=0,ncol=3)
    colnames(sub_0)<-c( "snp","phenotype","correlation")
    }
  return(sub_0)
  }
  },
  phen = phen,
  strict =strict,
  phen1 = phen1,
  strict1 = strict1,
  rel_genes = rel_genes,
  hap = hap,
  image = "eu.gcr.io/cytoreason/ci-cytoreason.single.cell-package:develop_0.2.0",
  memory_request='1Gi', 
  replace_image_tags = TRUE,
  force_execution = TRUE,
  simplify = FALSE)

temp <- get_outputs_dist("wf-882ac52e5d")
temp <- get_outputs_dist(res)


#wf-882ac52e5d- 100 with within gene control 0.85 
#wf-24d7a7fdfa- 100 with within gene control 0.84 (first step only) 

#wf-207f05f89d 100 with within gene control with replacement 0.75

#wf-6c34ad50ac - 1000 no gene control 0.78
#wf-3ce0a9f0d3 - 1000 with gene control 0.752

res_perm <- ldply(temp, function(lis){
  print(length(lis))
  if(length(lis)>0){
  if(colnames(lis$output.rds) %in% "snp"){
  ln <- length(unique(lis$output.rds$snp))}
  else {
    ln=0
  }
  ln
  }
})

res_perm <- ldply(temp, function(lis){
  print(length(lis))
  if(length(lis)>0){
    if(colnames(lis) %in% "snp"){
      ln <- length(unique(lis$snp))}
    else {
      ln=0
    }
    ln
  }
})

res_perm <-rbind(res_perm,data.frame(V1=rep(x = 0,9)))
colnames(res_perm)[2] <-"number_of_validated_snps"
ggplot(res_perm, aes(x=number_of_validated_snps))+geom_density()+scale_x_continuous(breaks=seq(0,3000,by=100))+ theme_bw()
ecdf(res_perm$number_of_validated_snps)(765)
ecdf(res_perm$V1)(1655)
