library("ggplot2")
install.packages("egg")
library("egg")

tot2 <- read.csv("capsule/code/CCanalysis/scripts/SignalPropagation/conservation_analysis_updated_adj_005.csv")
tot2 <- tot2 %>%
  select(Gene.name,  cons,phen, tr_expr, shared_united, ad_united)

conservation_avg01 <- read.csv("capsule/code/CCanalysis/scripts/Conservation/conservation_avg01.csv")
conservation_avg01 <- conservation_avg01 %>%
                    select(Gene.name, Mean)

conservation_avg01 <- conservation_avg01  %>%
                            dplyr::rename("cons"="Mean") 
conservation_avg01$"phen"<-"none"
conservation_avg01$"tr_expr"<-"none"
conservation_avg01$"shared_united"<-"none"
conservation_avg01$"ad_united"<-"none"
tot2 <- rbind(tot2,conservation_avg01)

tr <- ggplot(tot2, aes(cons, colour = as.factor(tr_expr))) +
  stat_ecdf(geom = "step")+ 
  theme_bw() + theme(text=element_text(size=15)) + 
  coord_cartesian(xlim = c(0, 350)) + labs(title="p.ks=8.176e-06")

ks.test(tot2$cons[which(tot2$tr_expr=="trans")],tot2$cons[which(tot2$tr_expr=="cis")],alternative = "greater")

cr <- ggplot(tot2, aes(cons, colour = as.factor(shared_united))) +
  stat_ecdf(geom = "step")+ theme_bw()+theme(text=element_text(size=15))+
  coord_cartesian(xlim = c(0, 350)) +labs(title="p.ks=1.064e-05")

ks.test(tot2$cons[which(tot2$shared_united=="shared")],tot2$cons[which(tot2$shared_united=="cell_specific")], alternative = "greater")

tot2_shared_united <- subset(tot2, shared_united %in% "shared")
pr <- ggplot(tot2_shared_united, aes(cons, colour = as.factor(ad_united))) +
  stat_ecdf(geom = "step")+ theme_bw()+theme(text=element_text(size=15)) +
  coord_cartesian(xlim = c(0, 350))+labs(title="p.ks=5.895e-05")

ks.test(tot2$cons[which(tot2$ad_united=="cis_shared")],tot2$cons[which(tot2$ad_united=="trans_shared")], alternative = "less")

tot2_cell_specific <- subset(tot2, shared_united %in% "cell_specific")
kr <- ggplot(tot2_cell_specific, aes(cons, colour = as.factor(ad_united))) +
  stat_ecdf(geom = "step")+ theme_bw()+theme(text=element_text(size=15)) +
  coord_cartesian(xlim = c(0, 350))+labs(title="p.ks=0.07828")

ks.test(tot2$cons[which(tot2$ad_united=="cis_cell_specific")],tot2$cons[which(tot2$ad_united=="trans_cell_specific")], alternative = "less")

figure <- ggarrange(tr, cr, pr,kr,
                    labels = c("a", "b", "c","d"),
                    ncol = 2, nrow = 2)
figure



