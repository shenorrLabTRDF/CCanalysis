
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
BiocManager::install("annotationTools")
BiocManager::install("VariantAnnotation")

BiocManager::install("fpc")
BiocManager::install("QTLRel")
BiocManager::install("regress")
BiocManager::install("RUnit")

library("BSgenome.Mmusculus.UCSC.mm10")
library("annotationTools")
library("VariantAnnotation")
library("fpc")
library("QTLRel")
library("regress")
library("RUnit")
load_all("capsule/code/DOQTL/")
library(abind)
library(plyr)
library(dplyr)

load("capsule/code/CCanalysis/scripts/mapping//MM_snps.Rdata")
load("capsule/code/CCanalysis/scripts/mapping/CCgenomesB38adj.Rdata")
outdir <- "capsule/code/CCanalysis//scripts/mapping/results/"
pheno <- read.csv("capsule/code/CCanalysis//scripts/mapping/phen_BM_updated_20200711.csv")

pheno <- pheno[order(pheno[,1]),]
ms_pheno <- c(11,38)
pheno <- pheno[-ms_pheno,]

intervals <- read.csv("capsule/code/CCanalysis//scripts/mapping/rel_intervals.csv")
intervals2 <- read.csv("capsule/code/CCanalysis/scripts/mapping//rel_intervals2.csv")
  
strict_matrix <- read.csv("capsule/code/CCanalysis/scripts/mapping/strict_matrix_old.csv")
strict_matrix <- strict_matrix %>%
                           dplyr::select(samples,Tcells,NKcells,ly6gPos,classicalMonocytes,proB,earlyB,lateB,CD8Tcells,CD4Tcells,StemCells,Bcells)

colnames(strict_matrix)[2:12] <- colnames(pheno)[3:13]

###here we are choosing i (3< i< 13),as it was an external variable indicating the column in pheno to anable job creation on the cluster
i=10 #"CD8Tcells"
pheno <- cbind(pheno[,1:2],pheno[,i,drop=FALSE])
pheno <- subset(pheno, Sample %in% strict_matrix$samples[which(strict_matrix[,i-1]!=0)])

model.probs <- model.probs[,,which(dimnames(model.probs)[3][[1]] %in% c(as.character(intervals$marker),intervals2$marker))]
m_sort <- abind(model.probs,model.probs,along=1)
m_sort <- m_sort[-which(rownames(m_sort) %in% c("CC016.GeniUnc", "CC021.Unc","CC057.Unc")),,] 
m_sort <- abind(m_sort,model.probs[which(rownames(model.probs) %in% c("CC016.GeniUnc", "CC021.Unc","CC057.Unc")),,],along=1)
#m_sort <- abind(m_sort,model.probs[which(rownames(model.probs) %in% c("CC009.Unc", "CC033.GeniUnc")),,],along=1)
m_sort <- m_sort[dimnames(m_sort)[[1]] %in% pheno$Sample,,]
m_sort <- m_sort[order(dimnames(m_sort)[1][[1]]),,]

MM_snps <- MM_snps[MM_snps$marker %in% c(as.character(intervals$marker),intervals2$marker),]

ref <- cbind(dimnames(m_sort)[1][[1]],paste0("F",seq(1:length(dimnames(m_sort)[1][[1]]))))
ref <- as.data.frame(ref)
dimnames(m_sort)[1][[1]]<-paste0("F",seq(1:length(dimnames(m_sort)[1][[1]])))

rownames(pheno) <- ref$V2
pheno[,1] <- rownames(pheno)
pheno[,1] <- as.factor(pheno[,1])
K <- kinship.probs(m_sort)


sub_pheno <- pheno
covar = data.frame(sex = as.numeric(sub_pheno$Sex == "M"))
rownames(covar) = rownames(sub_pheno)
qtl = scanone(pheno=sub_pheno, pheno.col = colnames(sub_pheno)[3], probs = m_sort,addcovar = covar, snps = MM_snps)
temp <- rbind(qtl$lod[[1]],qtl$lod[[2]])
temp <- subset(temp,temp$lod>2)

final <- matrix(ncol = 11, nrow = 0)
results <- do.call('rbind',lapply(1:nrow(pheno),function(rowType,strains,MM_snps,K,m_sort){
    strains <- pheno$Sample
    strains <- strains[-which(as.character(strains)==as.character((pheno[rowType,1])))]
    sub_pheno <- subset(pheno,as.character(pheno[,1]) %in% as.character(strains))
    covar = data.frame(sex = as.numeric(sub_pheno$Sex == "M"))
    rownames(covar) = rownames(sub_pheno)
    qtl = scanone(pheno=sub_pheno, pheno.col = colnames(sub_pheno)[3], probs = m_sort,addcovar = covar, snps = MM_snps)
    tempo <- rbind(qtl$lod[[1]],qtl$lod[[2]])
    temp2 <- cbind(tempo,colnames(pheno)[3],pheno[rowType,1])
    #final <- rbind(final,temp2[which(temp2$lod>2),])
    final <- rbind(final,temp2)
},strains,MM_snps,K,m_sort))


sub_x <- rbind(temp,results[,1:9])
len = nrow(pheno)+1

sub <- sub_x %>% dplyr::count(marker, sort = TRUE)

# results2 <- ddply(sub_x,.(marker),function(x){
# if (nrow(x) == len){
#       
#      x$m_lod <- mean(x$lod)
# 
#     }
# })
results2 <- subset(temp, marker %in% sub$marker[which(sub$n==len)])  

write.csv(results2,paste0(outdir,colnames(pheno)[3],".csv"))
