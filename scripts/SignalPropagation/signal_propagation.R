library("dplyr")
library("plyr")
library("factoextra")
library("pheatmap")
library("reshape2")
library("tidyr")


#validatedGenes <- read.csv("capsule/code/CCanalysis/scripts/validation/validated_genes.csv")
#genes <- read.csv("CC/scripts/validation/genes_with _direction_validation.csv")
genes <- read.csv("capsule/code/CCanalysis/scripts/validation/gene_cor_greater_0_updated.csv")
#genes <- read.csv("capsule/code/CCanalysis/scripts/validation/rel_genes_validated.csv")
genes <- genes %>%
           dplyr::rename("SNP"="marker")

all <- read.csv("capsule/code/CCanalysis/scripts/Figures/all2.csv")
imm_exons_intervals_rel <- read.csv("capsule/code/CCanalysis/scripts/Loci_filtering/rel_intervals_updated.csv")
all<- subset(all, marker %in% imm_exons_intervals_rel$marker)

RawAssociationMatrix <- read.csv("capsule/code/CCanalysis/scripts/SignalPropagation/all_associations.csv")
RawAssociationMatrix <- RawAssociationMatrix[-which(RawAssociationMatrix$phen %in% "preproB"),]
RawAssociationMatrix <- RawAssociationMatrix[which(RawAssociationMatrix$marker %in% (genes$SNP)),]


assM <- matrix(ncol= length(unique(RawAssociationMatrix$phen)), nrow=length(unique(RawAssociationMatrix$marker)))
rownames(assM) <- unique(RawAssociationMatrix$marker)
colnames(assM) <- unique(RawAssociationMatrix$phen)

assM <- apply(assM,c(1,2), function(x){x<-0})

for (r in 1:nrow(RawAssociationMatrix)){
  assM[RawAssociationMatrix$marker[r],RawAssociationMatrix$phen[r]]<-RawAssociationMatrix$lod[r]
  
}

assM <- as.data.frame(assM)

assN <- assM
for ( i in 1:nrow(assN)){
  for ( j in 1:9){
    if(assN[i,j] < 5){
      assN[i,j] <- 0
    }
  }
}
assN <- t(apply(assN,1,function(x){x/max(x)}))

# Silhouette method
res <- fviz_nbclust(assN, kmeans, method = "silhouette",k.max = 10)+
  labs(subtitle = "Silhouette method")
optNumC <- which(res$data$y %in% (max(res$data$y[1:10])))
ph <- pheatmap(assN,show_rownames = F,kmeans_k = 8,cluster_cols = F)

mem <- ph$kmeans$cluster
assN <- as.data.frame(assN)
assN$"mem" <- sapply(rownames(assN),function(x){mem[match(x,names(mem))]})
assN$"gene"<- rownames(assN)

assn_lf <- melt(assN,id.vars = c("mem","gene"))
assn_lf$mem <- as.character(assn_lf$mem)
ass_agg <-
  assn_lf %>% 
  group_by(mem, variable) %>% 
  dplyr::summarize(Mean = median(value, na.rm = TRUE)) %>% 
  group_by(mem) %>% 
  dplyr::mutate(pmedian = Mean / max(Mean))

#plot
ggplot(ass_agg, aes(x = pmedian)) + geom_density(alpha = 0.25) + scale_x_continuous(breaks = seq(0, 1, by =
                                                                                                   0.2)) + scale_y_continuous(breaks = seq(0, 3, by = 0.2)) + theme_linedraw() + theme_light()
cent <- ph$kmeans$centers
tot <- matrix(ncol=11,nrow=0)
for (c in 1:nrow(cent)){
  sub_tot <- rownames(assN)[which(assN$mem==c)]
  pops <- colnames(cent)[which(cent[c,] > 0.4)]
  for (pop in pops){
    sub_all <- all[which(all$phen==pop),] 
    sub_all <- sub_all[which(sub_all$marker %in% sub_tot),]  
   tot <- rbind(tot,sub_all)
  }
}

###adding genes taht couldn't be propagated
missed <-genes[which(!unique(genes$SNP) %in% unique(tot$marker)),] 
missed <- missed %>%
  dplyr::select(SNP,phen) %>%
  dplyr::rename("marker"="SNP")

##binding genes and propagated table
tot <- tot %>%
         dplyr::select(marker,phen)

tot <- rbind(tot,missed)

##adding gene names
#gene_intervals <- read.csv("capsule/code/CCanalysis/scripts/SignalPropagation/imm_exons_intervals.csv")
gene_intervals <- read.csv("capsule/code/CCanalysis/scripts/Loci_filtering/rel_intervals_updated.csv")

gene_intervals <- gene_intervals %>%
  dplyr::select(marker, Gene.name)

tot2 <- left_join(tot, gene_intervals,relationship="many-to-many")
tot2 <- tot2[!duplicated(tot2[,c("phen", "Gene.name","marker")]),] 

#sanity check, are there any genes that are more than in one locus
dat <- tot2 %>%                    # take the data.frame "data"
  dplyr::group_by(phen,Gene.name) %>%          # Then, with the filtered data, group it by "bb"
  dplyr::summarise(Unique_Elements = n_distinct(marker)) %>%   # Now summarise with unique elements per group
  ungroup()

#adding conservation scores
cons_scores <- read.csv("capsule/code/CCanalysis/scripts/SignalPropagation/conservationAllGenes.csv")
cons_scores <-subset(cons_scores, !is.na(avScore))
conservation_avg <-  cons_scores %>% 
  group_by(Gene.name) %>%
  dplyr::summarize(Mean = max(avScore, na.rm=TRUE))
tot2 <- left_join(tot2,conservation_avg)



###adding expression
expr <- read.csv("capsule/code/CCanalysis/scripts/Expression/expression_long_format.csv")
expr <- melt(expr)
colnames(expr) <- c("gene", "phen","value")
expr <- as.data.frame(expr)
expr$phen <- as.character(expr$phen)
expr$phen[which(expr$phen %in% "classicalMonocytes")] <- "Monocytes"
expr$phen[which(expr$phen %in% "ly6gPos")] <- "Granulocytes"
expr$phen[which(expr$phen %in% "restB")] <- "lateB"
expr$phen[which(expr$phen %in% "classicalMonocytes")] <- "Monocytes"

expr_B <- read.csv("capsule/code/CCanalysis//scripts/Expression/expression_long_format_Bcells.csv")
expr_B <- expr_B %>%
            dplyr::select(-X) %>%
          dplyr::rename("value"="Bcells")

expr <- rbind(expr,expr_B)
expr$phen[which(expr$phen %in% "StemCells.")] <- "StemCells"

tot2 <- left_join(tot2,expr, by= c("phen"="phen","Gene.name"="gene"))
tot2 <- tot2 %>%
  dplyr::rename("expr"="value") %>%
 dplyr::rename("cons"="Mean")


immgen_syn <- read.csv("capsule/code/CCanalysis/scripts/Loci_filtering/Filtered_variants/immgen_genes_synonyms.csv")
tot2_sub <- subset(tot2,is.na(tot2$expr))
tot2_sub$"syn" <- sapply(tot2_sub$Gene.name, function(x){immgen_syn$Gene.Synonym[match(x, immgen_syn$Gene.name)]})
tot2_sub$syn[which(tot2_sub$Gene.name %in% "Desi1")] <- "Pppde2"
tot2_sub$gene <- tot2_sub$syn
tot2_sub <- left_join(tot2_sub, expr)
tot2_sub <- tot2_sub[!is.na(tot2_sub$value),]
tot2_sub$expr <-tot2_sub$value
tot2_sub <- tot2_sub %>%
  dplyr::select(marker, phen, Gene.name,cons,expr)
tot2 <- tot2[!is.na(tot2$expr),]
tot2 <- rbind(tot2,tot2_sub)



tot2$expr[which(tot2$expr==0)]<-"trans"
tot2$expr[which(tot2$expr!="trans")]<-"cis"



tot2 <- tot2[!is.na(tot2$cons),]
tot2$"adj" <- paste0(tot2$Gene.name,"_",tot2$phen)
tot2 <- subset(tot2,!duplicated(tot2$adj))

tot2$gene<-tot2$Gene.name
tot2$"tr_expr"<-"cis"
tot2 <- ddply(tot2,.(gene), function(gen){
  if(any(unique(gen$expr) %in% "trans")){
    gen$tr_expr<-"trans"
  }
  gen
})

##adding information regarding association sharing
tot2$phen_united<-tot2$phen
tot2$phen_united[which(tot2$phen %in% c("lateB","proB"))]<-"Bcells"
tot2$"shared" <- "cell_specific"
tot2 <- ddply(tot2,.(gene), function(gen){
  if(length(unique(gen$phen))>1)
  {
    gen$shared <-"shared"
  }
  gen
})

tot2$"shared_united" <- "cell_specific"
tot2 <- ddply(tot2,.(gene), function(gen){
  if(length(unique(gen$phen_united))>1)
  {
    gen$shared_united <-"shared"
  }
  gen
})

tot2$ad <- paste0(tot2$tr_expr,"_",tot2$shared)
tot2_shared <- subset(tot2, shared %in% "shared")
ggplot(tot2_shared, aes(cons, colour = as.factor(ad))) +
  stat_ecdf(geom = "step")
tot2$ad_united <- paste0(tot2$tr_expr,"_",tot2$shared_united)

write.csv(tot2,"capsule/code/CCanalysis/scripts/SignalPropagation/conservation_analysis_updated_adj_005.csv")



