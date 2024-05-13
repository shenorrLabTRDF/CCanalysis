
library(ggplot2)
library("plyr")

validatedGenes <- read.csv("capsule/code/CCanalysis/scripts/SignalPropagation/conservation_analysis_updated_adj_005.csv")

validatedGenes_agg <- validatedGenes %>%
                          dplyr::group_by(phen,expr) %>%
                          dplyr::summarise(n=n_distinct(Gene.name))

cell_agg <- ddply(validatedGenes_agg, .(phen), function(phen){
  phen$n[which(phen$expr %in% "cis")]/phen$n[which(phen$expr %in% "trans")]
})

validatedGenes_agg$phen <- factor(validatedGenes_agg$phen, levels =  cell_agg$phen[order(cell_agg$V1, decreasing = TRUE)])

ggplot(validatedGenes_agg, aes(fill = expr, y=n, x=phen)) + 
  geom_bar(position="dodge", stat="identity") + theme_bw() +
  scale_y_continuous(breaks = seq(0,200,by=25))

write.csv(validatedGenes_agg,"validatedGenes_by_expr.csv")

###option with proportions
colnames(cell_agg)<-c("phenotype","trans")
cell_agg$"cis"<-100- cell_agg$trans
cell_agg_lf <- reshape2::melt(cell_agg, by = "phenotype")

cell_agg_lf$phenotype <- factor(cell_agg_lf$phenotype, levels =  cell_agg$phen[order(cell_agg$trans, decreasing = TRUE)])

ggplot(cell_agg_lf, aes(fill = variable, y=value, x=phenotype)) + 
  geom_bar(position="fill", stat="identity") + theme_bw()
