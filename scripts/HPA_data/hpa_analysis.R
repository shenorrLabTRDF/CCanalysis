library(tidyr)
library(dplyr)
library(plyr)
library("ggplot2")

##monocytes is a combination of all monocyte types
class_mono <- read.delim2("capsule/code/CCanalysis/scripts/HPA_data/expr_files/class_mono.tsv")
int_mono <- read.delim2("capsule/code/CCanalysis/scripts/HPA_data/expr_files/int_monocytes.tsv")
non_class_mono <- read.delim2("capsule/code/CCanalysis/scripts/HPA_data/expr_files/non_class_mono.tsv")
mono <- rbind(class_mono,int_mono, non_class_mono)

baso <- read.delim2("capsule/code/CCanalysis/scripts/HPA_data/expr_files/basophils_exp.tsv")
neu <- read.delim2("capsule/code/CCanalysis/scripts/HPA_data/expr_files/neutrophils_exp.tsv")
eos <- read.delim2("capsule/code/CCanalysis/scripts/HPA_data/expr_files/eosonophils_exp.tsv")

#lymphocytes is a combination of all B and T cell subsets
treg <- read.delim2("capsule/code/CCanalysis/scripts/HPA_data/expr_files/treg.tsv")
naivecd8 <- read.delim2("capsule/code/CCanalysis/scripts/HPA_data/expr_files/naiveCD8.tsv")
naiveCD4 <- read.delim2("capsule/code/CCanalysis/scripts/HPA_data/expr_files/naiveCD4.tsv")
memCD8 <- read.delim2("capsule/code/CCanalysis/scripts/HPA_data/expr_files/memoryCD8.tsv")
memCD4 <- read.delim2("capsule/code/CCanalysis/scripts/HPA_data/expr_files/memoryCD4.tsv")
naiveB <- read.delim2("capsule/code/CCanalysis/scripts/HPA_data/expr_files/naiveB.tsv")
memB <- read.delim2("capsule/code/CCanalysis/scripts/HPA_data/expr_files/memoryB.tsv")
gdT <- read.delim2("capsule/code/CCanalysis/scripts/HPA_data/expr_files/gdT.tsv")
maitT <- read.delim2("capsule/code/CCanalysis/scripts/HPA_data/expr_files/maitT.tsv")

lymphocytes <- rbind(treg, naivecd8, naiveCD4, memCD4, memCD8, naiveB, memB, gdT, maitT)

vuck_study <- read.csv("capsule/code/CCanalysis/scripts/Vuckovic_associated_genes.csv")
conservation <- read.csv("capsule/code/CCanalysis/scripts/HPA_data/vuck_study_conservation.csv")

vuck_study$"expr" <-"trans"

vuck_study <- vuck_study %>% 
  mutate(Ensembl.Gene.ID.s..for.Most.Serious.Consequence = strsplit(Ensembl.Gene.ID.s..for.Most.Serious.Consequence, ",")) %>% 
  unnest(Ensembl.Gene.ID.s..for.Most.Serious.Consequence)

vuck_study_mono <- subset(vuck_study, vuck_study$Associated.Blood.Index %in% "MONO%")
vuck_study_mono$expr[which(vuck_study_mono$Ensembl.Gene.ID.s..for.Most.Serious.Consequence %in% mono$Ensembl)] <- "cis"

vuck_study_eos<-subset(vuck_study, vuck_study$Associated.Blood.Index %in% "EO%")
vuck_study_eos$expr[which(vuck_study_eos$Ensembl.Gene.ID.s..for.Most.Serious.Consequence %in% eos$Ensembl)]<-"cis"

vuck_study_neu<-subset(vuck_study, vuck_study$Associated.Blood.Index %in% "NEUT%")
vuck_study_neu$expr[which(vuck_study_neu$Ensembl.Gene.ID.s..for.Most.Serious.Consequence %in% neu$Ensembl)]<-"cis"

vuck_study_baso<-subset(vuck_study, vuck_study$Associated.Blood.Index %in% "BASO%")
vuck_study_baso$expr[which(vuck_study_baso$Ensembl.Gene.ID.s..for.Most.Serious.Consequence %in% baso$Ensembl)]<-"cis"

vuck_study_lymph <-subset(vuck_study, vuck_study$Associated.Blood.Index %in% "LYMPH%")
vuck_study_lymph$expr[which(vuck_study_lymph$Ensembl.Gene.ID.s..for.Most.Serious.Consequence %in% lymphocytes$Ensembl)]<-"cis"


vuck_study <- rbind(vuck_study_baso,vuck_study_neu,vuck_study_eos,vuck_study_mono, vuck_study_lymph)

vuck_study$adj <- paste0(vuck_study$Associated.Blood.Index,"_",vuck_study$Ensembl.Gene.ID.s..for.Most.Serious.Consequence)
vuck_study <- vuck_study[-which(duplicated(vuck_study$adj)),]

##adding conservation information

#averaging on gene
conservation_avg <-  conservation %>% 
  group_by(Gene.stable.ID) %>%
  dplyr::summarize(Mean = max(avScore, na.rm=TRUE))

vuck_study$"conservation" <- sapply(vuck_study$Ensembl.Gene.ID.s..for.Most.Serious.Consequence,function(gene){conservation_avg$Mean[match(gene,conservation_avg$Gene.stable.ID)]})
vuck_study <- vuck_study[!is.na(vuck_study$conservation),]


vuck_study$"tr_expr"<-"cis"
vuck_study <- ddply(vuck_study,.(Ensembl.Gene.ID.s..for.Most.Serious.Consequence), function(gen){
  if(any(unique(gen$expr) %in% "trans")){
    gen$tr_expr<-"trans"
  }
  gen
})


##adding information regarding association sharing
vuck_study$"shared" <- "cell_specific"
vuck_study <- ddply(vuck_study,.(Ensembl.Gene.ID.s..for.Most.Serious.Consequence), function(gene){
  if(length(unique(gene$Associated.Blood.Index))>1)
  {
    gene$shared <-"shared"
  }
  gene
})

vuck_study$"lin"<-""
vuck_study$lin[which(vuck_study$Associated.Blood.Index %in% c("EO%","NEUT%","BASO%"))]<-"granulocytes"
vuck_study$lin[which(vuck_study$Associated.Blood.Index %in% c("LYMPH%"))]<-"lymphocytes"
vuck_study$lin[which(vuck_study$Associated.Blood.Index %in% c("MONO%"))]<-"monocytes"

vuck_study$"shared_united" <- "cell_specific"
vuck_study <- ddply(vuck_study,.(Ensembl.Gene.ID.s..for.Most.Serious.Consequence), function(gene){
  if(length(unique(gene$lin))>1)
  {
    gene$shared_united <-"shared"
  }
  gene
})

vuck_study$"property1" <- paste0(vuck_study$shared,"_",vuck_study$tr_expr)
vuck_study$"property2" <- paste0(vuck_study$shared_united,"_",vuck_study$tr_expr)

write.csv(vuck_study,"capsule/code/CCanalysis/conservation_analysis_vuck_study.csv")



###Chen study
chen_study <- read.csv("capsule/code/CCanalysis/scripts/Chen_association_data.csv")
chen_study <- chen_study %>% 
  mutate(SYMBOL = strsplit(SYMBOL, ",")) %>% 
  unnest(SYMBOL)
chen_study <- chen_study[!is.na(chen_study$SYMBOL),]

chen_study$"expr" <-"trans"
 

chen_study_mono <- subset(chen_study, chen_study$Trait %in% "MONO")
chen_study_mono$expr[which(chen_study_mono$SYMBOL %in% mono$Gene)]<-"cis"

chen_study_eos <- subset(chen_study, chen_study$Trait %in% "EOS")
chen_study_eos$expr[which(chen_study_eos$SYMBOL %in% eos$Gene)]<-"cis"

chen_study_neu <- subset(chen_study, chen_study$Trait %in% "NEU")
chen_study_neu$expr[which(chen_study_neu$SYMBOL %in% neu$Gene)]<-"cis"

chen_study_baso <- subset(chen_study, chen_study$Trait %in% "BASO")
chen_study_baso$expr[which(chen_study_baso$SYMBOL %in% baso$Gene)]<-"cis"

chen_study_lymph <-subset(chen_study, chen_study$Trait %in% "LYM")
chen_study_lymph$expr[which(chen_study_lymph$SYMBOL %in% lymphocytes$Gene)]<-"cis"


chen_study <- rbind(chen_study_baso,chen_study_neu,chen_study_eos,chen_study_mono, chen_study_lymph)

chen_study$adj <- paste0(chen_study$Trait,"_",chen_study$SYMBOL)
chen_study <- chen_study[-which(duplicated(chen_study$adj)),]


##adding conservation information
conservation <- read.csv("capsule/code/CCanalysis/scripts/HPA_data/chen_study_conservation.csv")

#averaging on gene
conservation_avg <-  conservation %>% 
  group_by(Gene.name) %>%
  dplyr::summarize(Mean = max(avScore, na.rm=TRUE))

chen_study$"conservation" <- sapply(chen_study$SYMBOL,function(gene){conservation_avg$Mean[match(gene,conservation_avg$Gene.name)]})
chen_study$conservation[which(is.na(chen_study$conservation))] <- sapply(chen_study$SYMBOL[which(is.na(chen_study$conservation))], function(x){vuck_study$conservation[match(x,vuck_study$Gene.Symbol.s..for.Most.Serious.Consequence)]})
chen_study <- chen_study[!is.na(chen_study$conservation),]


##adding information regarding association sharing

chen_study$"tr_expr"<-"cis"
chen_study <- ddply(chen_study,.(SYMBOL), function(gen){
  if(any(unique(gen$expr) %in% "trans")){
    gen$tr_expr<-"trans"
  }
  gen
})


chen_study$"shared" <- "cell_specific"
chen_study <- ddply(chen_study,.(SYMBOL), function(gene){
  if(length(unique(gene$Trait))>1)
  {
    gene$shared <-"shared"
  }
  gene
})

chen_study$"lin"<-""
chen_study$lin[which(chen_study$Trait %in% c("EOS","NEU","BASO"))]<-"granulocytes"
chen_study$lin[which(chen_study$Trait %in% c("LYM"))]<-"lymphocytes"
chen_study$lin[which(chen_study$Trait %in% c("MONO"))]<-"monocytes"

chen_study$"shared_united" <- "cell_specific"
chen_study <- ddply(chen_study,.(SYMBOL), function(gene){
  if(length(unique(gene$lin))>1)
  {
    gene$shared_united <-"shared"
  }
  gene
})

chen_study$"property1" <- paste0(chen_study$shared,"_",chen_study$tr_expr)
chen_study$"property2" <- paste0(chen_study$shared_united,"_",chen_study$tr_expr)

write.csv(chen_study,"capsule/code/CCanalysis/conservation_analysis_chen_study.csv")

###combined Vuck and Chen studies

vuck_filt <- vuck_study %>% select(Associated.Blood.Index, Gene.Symbol.s..for.Most.Serious.Consequence, conservation, expr,tr_expr,
                                   shared,lin,shared_united,property1,property2)
colnames(vuck_filt)[1:2] <- c("trait","symbol")

vuck_filt <- vuck_filt %>% mutate(trait = stringr::str_replace(trait, "EO%", "EOS"))
vuck_filt <- vuck_filt %>% mutate(trait = stringr::str_replace(trait, "NEUT%", "NEU"))
vuck_filt <- vuck_filt %>% mutate(trait = stringr::str_replace(trait, "LYMPH%", "LYM"))
vuck_filt <- vuck_filt %>% mutate(trait = stringr::str_replace(trait, "BASO%", "BASO"))
vuck_filt <- vuck_filt %>% mutate(trait = stringr::str_replace(trait, "MONO%", "MONO"))
vuck_filt$study <-"vuck"

chen_filt <- chen_study %>% select(Trait, SYMBOL, conservation, expr,tr_expr,
                                   shared,lin,shared_united,property1,property2)
colnames(chen_filt)[1:2]<-c("trait","symbol")
chen_filt$study <-"chen"

combined <- rbind(vuck_filt,chen_filt)
combined$adj <- paste0(combined$trait,"_",combined$symbol)
combined <- combined[-which(duplicated(combined$adj)),]

write.csv(combined,"capsule/code/CCanalysis/conservation_analysis_combined.csv")



###Orru study

orru_study <- read.csv("CC/scripts/orru_study_associations.csv")
conservation <- read.csv("CC/scripts/HPA_data/orru_study_conservation.csv")

orru_study$"expr" <- 0
#vuck_study<-vuck_study[-which(vuck_study$Ensembl.Gene.ID.s..for.Most.Serious.Consequence %in% ""),]

orru_study <- orru_study %>% 
  mutate(Gene.annotation....1Mb. = strsplit(Gene.annotation....1Mb., ",")) %>% 
  unnest(Gene.annotation....1Mb.)


##adding conservation information

conservation_avg <-  conservation %>% 
  group_by(Gene.name) %>%
  dplyr::summarize(Mean = max(avScore, na.rm=TRUE))

orru_study$"conservation" <- sapply(orru_study$Gene.annotation....1Mb.,function(gene){conservation_avg$Mean[match(gene,conservation_avg$Gene.name)]})
orru_study<-orru_study[!is.na(orru_study$conservation),]

##adding information regarding association sharing
orru_study$"shared" <- "cell_specific"
orru_study <- ddply(orru_study,.(Gene.annotation....1Mb.), function(gene){
  if(length(unique(gene$Trait))>1)
  {
    gene$shared <-"shared"
  }
  gene
})

ggplot(orru_study, aes(conservation, colour = shared)) +
  stat_ecdf(geom = "step")



chen_study_mono <- subset(chen_study, chen_study$Trait %in% "MONO")
chen_study_mono$expr[which(chen_study_mono$SYMBOL %in% mono$Gene)]<-1

chen_study_eos <- subset(chen_study, chen_study$Trait %in% "EOS")
chen_study_eos$expr[which(chen_study_eos$SYMBOL %in% eos$Gene)]<-1

chen_study_neu <- subset(chen_study, chen_study$Trait %in% "NEU")
chen_study_neu$expr[which(chen_study_neu$SYMBOL %in% neu$Gene)]<-1

chen_study_baso <- subset(chen_study, chen_study$Trait %in% "BASO")
chen_study_baso$expr[which(chen_study_baso$SYMBOL %in% baso$Gene)]<-1

chen_study_lymph <-subset(chen_study, chen_study$Trait %in% "LYM")
chen_study_lymph$expr[which(chen_study_lymph$SYMBOL %in% lymphocytes$Gene)]<-1


chen_study <- rbind(chen_study_baso,chen_study_neu,chen_study_eos,chen_study_mono, chen_study_lymph)

