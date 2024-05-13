

all_associations <- read.csv("capsule/code/CCanalysis/scripts/Figures/all2.csv")
all_associations <- all_associations %>%
  dplyr::filter(!phen %in% "preproB")

#fdr values taken from previous aqnalysis
fdr_key <- data.frame(cell_types = c("Bcells","CD4Tcells","CD8Tcells","Granulocytes","lateB","Monocytes","NKcells","proB","StemCells"),
                      fdr_tr = c(5.75,7,8.5,6,6,6,8.9,6.6,10.11)) 


all_associations <- ddply(all_associations,.(phen),function(phenotype){
                       phenotype <- subset(phenotype, lod > fdr_key$fdr_tr[fdr_key$cell_types %in% phenotype$phen])
                       phenotype
  
})

imm_exons_intervals_rel <- read.csv("capsule/code/CCanalysis/scripts/Loci_filtering/rel_intervals_updated.csv")
imm_exons_intervals_rel$"adj"<- paste0(imm_exons_intervals_rel$marker,"_",imm_exons_intervals_rel$Gene.name)
imm_exons_intervals_rel <- imm_exons_intervals_rel[which(!duplicated(imm_exons_intervals_rel$adj)),]
imm_exons_intervals_rel <- imm_exons_intervals_rel %>%
                           dplyr::filter(marker %in% all_associations$marker) %>%
                           dplyr::select(marker,Gene.name)
all_associations <-left_join(imm_exons_intervals_rel,all_associations, by = "marker", relationship="many-to-many")

write.csv(all_associations , "capsule/code/CCanalysis/scripts/Loci_filtering/rel_genes_remapped.csv")
