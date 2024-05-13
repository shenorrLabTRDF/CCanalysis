library("dplyr")
library("plyr")

gene_assignment <- read_excel("capsule/code/CCanalysis/scripts/Loci_filtering/Filtered_variants/gene_assignment.xls")
immgen_genes_synonyms <- read.delim2("capsule/code/CCanalysis/scripts/Loci_filtering/Filtered_variants/MGIBatchReport_immgen.txt")
SNPs_outfiltered_tot <- read.csv("capsule/code/CCanalysis/scripts/Loci_filtering/Filtered_variants/SNPs_outfiltered_tot.csv")
Indels_outfiltered_tot <- read.csv("capsule/code/CCanalysis/scripts/Loci_filtering/Filtered_variants/Indels_outfiltered_tot.csv")
SVs_outfiltered_tot <- read.csv("capsule/code/CCanalysis/scripts/Loci_filtering/Filtered_variants/SVs_outfiltered_tot.csv")

immgen_genes_synonyms_old<- subset(immgen_genes_synonyms, Input.Type %in% "old symbol")
imm_exons_intervals <- read.csv("capsule/code/CCanalysis/scripts/Loci_filtering/Filtered_variants/imm_exons_intervals.csv")
imm_exons_intervals <- subset(imm_exons_intervals, !duplicated(Gene.name))
gene_assignment$"intervals" <- 0
gene_assignment$intervals[which(tolower(gene_assignment$Gene) %in% tolower(imm_exons_intervals$Gene.name))]<-1
gene_assignment$"old"<- sapply(gene_assignment$Gene, function(x){immgen_genes_synonyms_old$Symbol[match(x,immgen_genes_synonyms_old$Input)]})
gene_assignment$"intervals2" <- 0
gene_assignment$intervals2[which(tolower(gene_assignment$old) %in% tolower(imm_exons_intervals$Gene.name))]<-1
gene_assignment <- gene_assignment %>%
                       dplyr::rowwise() %>%
                       dplyr::mutate(SUM =sum(c(intervals,intervals2)))

gene_assignment$Gene[which(gene_assignment$SUM %in% 1 & gene_assignment$intervals2 %in% 1)] <- gene_assignment$old[which(gene_assignment$SUM %in% 1 & gene_assignment$intervals2 %in% 1)]

imm_exons_intervals$"immgen_gene" <- 0
imm_exons_intervals$immgen_gene[which(tolower(imm_exons_intervals$Gene.name) %in% tolower(gene_assignment$Gene))]<-1
imm_exons_intervals_rel <- subset(imm_exons_intervals,immgen_gene %in% 1)

imm_exons_intervals$"snps" <- 0
imm_exons_intervals$snps[which(tolower(imm_exons_intervals$Gene.name) %in% tolower(SNPs_outfiltered_tot$Gene.name))]<-1

imm_exons_intervals$"indels" <- 0
imm_exons_intervals$indels[which(tolower(imm_exons_intervals$Gene.name) %in% tolower(Indels_outfiltered_tot$Gene.name))]<-1

imm_exons_intervals$"svs" <- 0
imm_exons_intervals$svs[which(tolower(imm_exons_intervals$Gene.name) %in% tolower(SVs_outfiltered_tot$gene))]<-1

imm_exons_intervals_full <- read.csv("capsule/code/CCanalysis/scripts/Loci_filtering/Filtered_variants/imm_exons_intervals.csv")

imm_exons_intervals_rel <- subset(imm_exons_intervals,immgen_gene %in% 1)


imm_exons_intervals_rel <- imm_exons_intervals_rel %>%
  dplyr::rowwise() %>%
  dplyr::mutate(total= sum(c(snps,indels,svs)))

imm_exons_intervals_rel <- imm_exons_intervals_rel %>%
                           dplyr::filter( !total %in% 0)
imm_exons_intervals_rel <- subset(imm_exons_intervals_full, Gene.name %in% imm_exons_intervals_rel$Gene.name) 
write.csv(imm_exons_intervals_rel,"capsule/code/CCanalysis/scripts/Loci_filtering/rel_intervals_updated.csv")


