install.packages("lattice")
library(lattice)
library("dplyr")
library("ggplot2")
install.packages("ggtext")
library("ggtext")
source("~/capsule/code/CCanalysis/scripts/Figures/manhattan_function.R")

all_associations <- read.csv("capsule/code/CCanalysis/scripts/Figures/all2.csv")
rel_genes <- read.csv("capsule/code/CCanalysis/scripts/Figures/rel_genes_ALL_FDRfiltered2_beforeImp.csv")


imm_exons_intervals <- read.csv("capsule/code/CCanalysis/scripts/SignalPropagation/imm_exons_intervals.csv")
imm_exons_intervals$"adj"<- paste0(imm_exons_intervals$marker,"_",imm_exons_intervals$Gene.name)
imm_exons_intervals <- imm_exons_intervals[which(!duplicated(imm_exons_intervals$adj)),]

###expand rel_genes using all the intervals for a specific gene
imm_exons_intervals_rel <- subset(imm_exons_intervals,Gene.name %in% rel_genes$Marker.Symbol)
marker_gene_pairs <- imm_exons_intervals_rel %>% dplyr::select(marker,Gene.name)
rel_genes <-left_join(marker_gene_pairs, rel_genes,by=c("Gene.name"="Marker.Symbol"),relationship="many-to-many")

rel_genes$"adj" <- paste0(rel_genes$marker,"_",rel_genes$phen)

all_associations <- all_associations %>%
                      dplyr::filter(!phen %in% "preproB")
mapfile <- read.csv("capsule/code/CCanalysis/scripts/Figures/mapfile.csv")
conservation_analysis <- read.csv("capsule/code/CCanalysis/scripts/validation/gene_cor_greater_0.csv")
imm_exons_intervals_conseravation_analysis <-  subset(imm_exons_intervals,Gene.name %in% conservation_analysis$Marker.Symbol)
mg_pairs <- imm_exons_intervals_conseravation_analysis %>%
                            dplyr::select(marker,Gene.name)
conservation_analysis <-left_join(mg_pairs,conservation_analysis, by = c("Gene.name"="Marker.Symbol"), relationship="many-to-many")

conservation_analysis$"adj" <- paste0(conservation_analysis$SNP,"_",conservation_analysis$phen)
all_associations <- left_join(all_associations,mapfile)
all_associations$"sig" <- 0
all_associations$"adj" <- paste0(all_associations$marker,"_",all_associations$phen)


all_associations$sig[which(all_associations$adj %in% conservation_analysis$adj)] <-1
# all_associations_max <- all_associations %>%
#                         dplyr::filter(sig==0) %>%
#                         dplyr::group_by(marker) %>%
#                         dplyr::filter(lod==max(lod))
# all_associations_bind <- rbind(all_associations[which(all_associations$sig==1),],all_associations_max)

all_associations$p[which(all_associations$sig==1)]<-sapply(all_associations$adj[which(all_associations$sig==1)], function(x){rel_genes$p[match(x, rel_genes$adj)]})
all_associations$sig[which(all_associations$sig==1)]<-all_associations$phen[which(all_associations$sig==1)]

ass_to_remove <- rel_genes$adj[which(!rel_genes$adj %in% conservation_analysis$adj)]
all_associations <- all_associations[-which(all_associations$adj %in% ass_to_remove),] 
#all_associations <- all_associations[sort(all_associations$chr,decreasing = FALSE),]
rm <- which(all_associations$lod>10 & all_associations$sig==0)
exp <- all_associations[rm,]
exp <- left_join(exp,imm_exons_intervals,relationship = "many-to-many", by= "marker")


ann <- all_associations$sig

pdf("tty_man3.pdf")

manhattan.plot(all_associations$chromosome, all_associations$bp, all_associations$p, annotate=list(ann,"Bcells"=list(col="#221ED6"),"lateB"=list(col="#4062AD"),
                                                                       "StemCells"=list(col="#AAA6D2"),"Monocytes"=list(col="#F04C5A"),"CD4Tcells"=list(col="#D09FC8"),
                                                                       "Granulocytes"=list(col="#E47825"),"0"=list(col="#d3d3d3")))
dev.off()

###updated version
all_associations <- read.csv("capsule/code/CCanalysis/scripts/Figures/all2.csv")
all_associations <- all_associations %>%
  dplyr::filter(!phen %in% "preproB")
imm_exons_intervals_rel <- read.csv("capsule/code/CCanalysis/scripts/Loci_filtering/rel_intervals_updated.csv")
all_associations<- subset(all_associations, marker %in% imm_exons_intervals_rel$marker)


mapfile <- read.csv("capsule/code/CCanalysis/scripts/Figures/mapfile.csv")
all_associations <- left_join(all_associations,mapfile)
 #rel_genes2 <- read.csv("capsule/code/CCanalysis/scripts/SignalPropagation/conservation_analysis_updated_adj_005.csv")
 rel_genes2 <- read.csv("capsule/code/CCanalysis/scripts/validation/gene_cor_greater_0_updated.csv")

 rel_genes2$"adj" <- paste0(rel_genes2$marker,"_",rel_genes2$phen)
all_associations$"adj" <- paste0(all_associations$marker,"_",all_associations$phen)


rel_genes <- read.csv("capsule/code/CCanalysis/scripts/validation/rel_genes_validated_updated.csv")
rel_genes$"adj" <- paste0(rel_genes$marker,"_",rel_genes$phen)

ass_to_remove <- rel_genes$adj[which(!rel_genes$adj %in% rel_genes2$adj)]
all_associations <- subset(all_associations,!(all_associations$adj %in% ass_to_remove)) 

all_associations$"sig" <- 0
all_associations$sig[which(all_associations$adj %in% rel_genes2$adj)] <-1
all_associations$sig[which(all_associations$sig==1)]<-all_associations$phen[which(all_associations$sig==1)]
all_associations$chr[which(all_associations$chr %in% "X")] <- "20"
all_associations$chr <- as.numeric(all_associations$chr) 
all_associations <- all_associations[order(all_associations$chr,decreasing = FALSE),]

ann <- all_associations$sig

pdf("tty_man3.pdf")

manhattan.plot(all_associations$chr, all_associations$bp, all_associations$p, annotate=list(ann,"Bcells"=list(col="#60AE47"),"lateB"=list(col="#4AAD95"),
                                                                                                   "StemCells"=list(col="#D1964D"),"Monocytes"=list(col="#EA96A3"),"CD4Tcells"=list(col="#999933"),
                                                                                                   "Granulocytes"=list(col="#4DABB3"),"proB"=list(col="#6DAEE2"),"CD8Tcells"=list(col="#C4A1EA"),
                                                                                                    "0"=list(col="#d3d3d3")))
dev.off()
