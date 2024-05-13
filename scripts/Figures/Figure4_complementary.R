library("dplyr")
library("plyr")
library(readxl)
install.packages("egg")
library("egg")

cons_scores <- read.csv("capsule/code/CCanalysis/scripts/SignalPropagation/conservationAllGenes.csv")
#cons_scores <- read.csv("capsule/code/CCanalysis/scripts/Conservation/mouse_conservation_scores_all.csv")

cons_scores <-subset(cons_scores, !is.na(avScore))
conservation_avg <-  cons_scores %>% 
  group_by(Gene.name) %>%
  dplyr::summarize(Mean = max(avScore, na.rm=TRUE))

gene_assignment <- read_excel("capsule/code/CCanalysis/scripts/Loci_filtering/Filtered_variants/gene_assignment.xls")
imm_exons_intervals_rel <- read.csv("capsule/code/CCanalysis/scripts/Loci_filtering/rel_intervals_updated.csv")
imm_exons_intervals_rel <- subset(imm_exons_intervals_rel, !duplicated(Gene.name))
immgen_genes_synonyms <- read.delim2("capsule/code/CCanalysis/scripts/Loci_filtering/Filtered_variants/MGIBatchReport_immgen.txt")
immgen_genes_synonyms_old<- subset(immgen_genes_synonyms, Input.Type %in% "old symbol")


gene_assignment$"intervals" <- 0
gene_assignment$intervals[which(tolower(gene_assignment$Gene) %in% tolower(imm_exons_intervals_rel$Gene.name))]<-1
gene_assignment$"old"<- sapply(gene_assignment$Gene, function(x){immgen_genes_synonyms_old$Symbol[match(x,immgen_genes_synonyms_old$Input)]})
gene_assignment$"intervals2" <- 0
gene_assignment$intervals2[which(tolower(gene_assignment$old) %in% tolower(imm_exons_intervals_rel$Gene.name))]<-1
gene_assignment <- gene_assignment %>%
  dplyr::rowwise() %>%
  dplyr::mutate(SUM =sum(c(intervals,intervals2)))

gene_assignment$Gene[which(gene_assignment$SUM %in% 1 & gene_assignment$intervals2 %in% 1)] <- gene_assignment$old[which(gene_assignment$SUM %in% 1 & gene_assignment$intervals2 %in% 1)]

gene_assignment$"in_set" <- 0
gene_assignment$in_set[which(tolower(gene_assignment$Gene) %in% tolower(imm_exons_intervals_rel$Gene.name))] <-1
gene_assignment$Gene <- tolower(gene_assignment$Gene)
conservation_avg$Gene.name <- tolower(conservation_avg$Gene.name)

conservation_avg <- left_join(conservation_avg,gene_assignment, by=c("Gene.name"="Gene"),relationship="many-to-many")
conservation_avg$in_set[which(is.na(conservation_avg$in_set))] <- 2
conservation_avg01 <- subset(conservation_avg, in_set %in% c(0,1))
conservation_avg12 <- subset(conservation_avg, in_set %in% c(1,2))



tr <- ggplot(conservation_avg, aes(Mean, colour = as.factor(in_set))) +
  stat_ecdf(geom = "step")+ theme_bw()+theme(text=element_text(size=15)) +
  coord_cartesian(xlim = c(0, 350))


pr <- ggplot(conservation_avg01, aes(Mean, colour = as.factor(in_set))) +
  stat_ecdf(geom = "step")+ theme_bw()+theme(text=element_text(size=15)) +
 labs(title="p.ks = 2.2e-16")
ks.test(conservation_avg01$Mean[which(conservation_avg01$in_set==0)],conservation_avg01$Mean[which(conservation_avg01$in_set==1)])


kr <- ggplot(conservation_avg12, aes(Mean, colour = as.factor(in_set))) +
  stat_ecdf(geom = "step")+ theme_bw()+theme(text=element_text(size=15)) 
ks.test(conservation_avg12$Mean[which(conservation_avg12$in_set==1)],conservation_avg12$Mean[which(conservation_avg12$in_set==2)])

figure <- ggarrange( tr,pr,kr,
                    labels = c("a", "b","c"),
                    ncol = 2, nrow = 2)
figure

write.csv(conservation_avg01,"capsule/code/CCanalysis/scripts/Conservation/conservation_avg01.csv")
conservation_avg$conservation_ex[which(!(conservation_avg$conservation ==0))]<-0