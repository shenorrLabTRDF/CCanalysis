####circular heatmap figure 4b
library(reshape2)
library(ggplot2)
library("pheatmap")

validatedGenes <- read.csv("capsule/code/CCanalysis/scripts/SignalPropagation/conservation_analysis_updated_adj_005.csv")
validatedGenes$"adj" <- paste0(validatedGenes$phen,"_",validatedGenes$Gene.name)
nData <- matrix(ncol = 3,nrow = length(unique(validatedGenes$Gene.name))*length(unique(validatedGenes$phen)))

colnames(nData) <- c('sample','gene','relevant')
nData <- as.data.frame(nData)
nData$sample <- unique(validatedGenes$phen)
nData$gene <- unique(validatedGenes$Gene.name)

addRows <- matrix(ncol=3,nrow=8*25)
colnames(addRows) <- colnames(nData)
nData <- rbind(nData,addRows)
nData$gene[2169:2368] <- sample(25,replace = F)

nData$sample[2169:2368] <- unique(validatedGenes$phen)
nData <- nData[order(nData$sample),]
nData$relevant <- 0
nData$"adj" <- paste0(nData$sample,"_",nData$gene)

###adding information about relevant gene-phen pairs
nData$relevant[which(nData$adj %in% validatedGenes$adj)]<-1
#setting trans as a different level
nData$relevant[which(nData$adj %in% validatedGenes$adj[which(validatedGenes$expr %in% "trans")])]<-2

newM <- dcast(nData,gene~sample,value.var = "relevant")
rownames(newM) <- newM[,1]
newM <- newM[,-1]

newM <- newM[,c("CD4Tcells","CD8Tcells","Granulocytes","Monocytes","Bcells","lateB","proB","StemCells")]


ph=pheatmap(newM,show_rownames = F,cluster_cols = FALSE)
nData<-as.data.frame(nData)
nData$gene<-as.factor(nData$gene)
nData$gene<-factor(nData$gene,levels=(rownames(newM)[ph$tree_row$order]))
nData$sample<-factor(nData$sample,levels=c(colnames(newM)))
nData$var2 = as.numeric(nData$sample) + 15
nData$var2<-nData$var2-10

svg('withModules2.svg')
ggplot(nData, aes(x=gene, y=var2,fill=relevant)) +
  geom_tile(colour="grey") +
  ylim(c(0, 50)) +
  coord_polar(theta="x") +
  theme(panel.background=element_blank(),
        axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_text(size=5))+ theme(axis.title.y=element_blank(),
                                                 axis.text.y=element_blank(),
                                                 axis.ticks.y=element_blank())
dev.off()


###plotting the version with cis/trans


nData$"expr" <- sapply(nData$adj, function(x){validatedGenes$expr[match(x,validatedGenes$adj)]})
nData$expr[which(is.na(nData$expr))] <- "none"
nData$"forPlot"<-paste0(nData$relevant,"_",nData$expr)
ggplot(nData, aes(x=gene, y=var2,fill=forPlot)) +
  geom_tile(colour="grey") +
  ylim(c(0, 50)) +
  coord_polar(theta="x") +
  theme(panel.background=element_blank(),
        axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_text(size=5))+ theme(axis.title.y=element_blank(),
                                                 axis.text.y=element_blank(),
                                                 axis.ticks.y=element_blank())
