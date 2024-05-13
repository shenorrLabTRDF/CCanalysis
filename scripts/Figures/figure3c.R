###associated genes
library("dplyr")
library("plyr")
library("reshape2")

genes <- read.csv("capsule/code/CCanalysis/scripts/SignalPropagation/conservation_analysis_updated_adj_005.csv")


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

##removing non expressed genes
expr_matrix<- expr %>%
     dplyr::filter(!value==0)
  
  
expr_matrix_agg <-  expr_matrix %>%
                   group_by(gene) %>%
  

## sub_imm gene expression
all_phen <- expand.grid(unique(genes$phen), unique(genes$phen))
colnames(all_phen) <- c("ass","exp")
tot <- matrix(ncol=8 ,nrow = 8)
rownames(tot) <- unique(genes$phen)
colnames(tot) <- unique(genes$phen)

tot <- apply(tot, c(1,2),function(x){x<-0})
phens <- unique(genes$phen)

for (row in 1:length(phens))
{
  assc <- genes$gene[genes$phen==phens[row]]
  expr <- genes$gene[genes$phen==phens[row] & genes$expr %in% "cis"]
  tot[as.character(phens[row]),as.character(phens[row])] <- length(which(as.character(assc) %in% expr))#/length(assc)
  trans <- assc[which(!(as.character(assc) %in% expr))]
  
  for(phen in phens){
    if(as.character(phen)!=as.character(phens[row])){
      expr2 <- expr_matrix$gene[which(expr_matrix$phen %in% phen)]
      
      tot[as.character(phens[row]),as.character(phen)] <- length(which(as.character(trans) %in% expr2))#/length(assc)
    }
  }
}

tot<- tot*100
tot <- round(tot,digits = 0)
write.csv(tot,"capsule/code/CCanalysis/circos_updated_absolute_values.csv")

diag(tot) <-0
write.csv(tot,"capsule/code/CCanalysis/circos_updated_absolute_values_trans.csv")
