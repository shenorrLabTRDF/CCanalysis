library("ggplot2")
install.packages("egg")
library("egg")


tot2 <- read.csv("capsule/code/CCanalysis/conservation_analysis_combined.csv")


###plot the data
tr <- ggplot(tot2, aes(conservation, colour = as.factor(tr_expr))) +
  stat_ecdf(geom = "step") + coord_cartesian(xlim = c(0, 350))+
  theme_bw() + theme(text=element_text(size=15)) + labs(title="5.992e-10")

ks.test(tot2$conservation[which(tot2$tr_expr=="trans")],tot2$conservation[which(tot2$tr_expr=="cis")],alternative = "greater")


cr <- ggplot(tot2, aes(conservation, colour = as.factor(shared))) +
  stat_ecdf(geom = "step")+ theme_bw()+theme(text=element_text(size=15))+coord_cartesian(xlim = c(0, 350))+
  labs(title="p.ks=0.5391")

ks.test(tot2$conservation[which(tot2$shared=="shared")],tot2$conservation[which(tot2$shared=="cell_specific")], alternative = "greater")


tot2_shared <- subset(tot2, shared %in% "shared")
pr <- ggplot(tot2_shared, aes(conservation, colour = as.factor(property1))) +coord_cartesian(xlim = c(0, 350))+
  stat_ecdf(geom = "step")+ theme_bw()+theme(text=element_text(size=15)) +
  labs(title="p.ks=0.0004983")

ks.test(tot2$conservation[which(tot2$property1=="shared_cis")],tot2$conservation[which(tot2$property1=="shared_trans")], alternative = "less")

tot2_cell_specific <- subset(tot2, shared %in% "cell_specific")
kr <- ggplot(tot2_cell_specific, aes(conservation, colour = as.factor(property1))) +
  stat_ecdf(geom = "step")+ theme_bw()+theme(text=element_text(size=15)) +
  coord_cartesian(xlim = c(0, 350))+labs(title="p.ks=1.286e-07")

ks.test(tot2$conservation[which(tot2$property1=="cell_specific_cis")],tot2$conservation[which(tot2$property1=="cell_specific_trans")], alternative = "less")

figure <- ggarrange(tr, cr, pr,kr,
                    labels = c("a", "b", "c","d"),
                    ncol = 2, nrow = 2)
figure

