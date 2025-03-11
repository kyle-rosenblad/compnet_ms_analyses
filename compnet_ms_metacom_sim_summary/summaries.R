options(scipen=999)
library(ggplot2)
library(ggthemes)
library(patchwork)

uc <- readRDS("summary_uncorr.RDS")
ucpa <- readRDS("summary_uncorr_presabs.RDS")
corr <- readRDS("summary_corr.RDS")
corrpa <- readRDS("summary_corr_presabs.RDS")

tableau_colors <- ggthemes::tableau_color_pal("Tableau 10")(10)  # Adjust number as needed

nullmod_colors <- c("compnet" = tableau_colors[6], 
               "c1" = tableau_colors[1], 
               "c2" = tableau_colors[2],
               "c3" = tableau_colors[3], 
               "c4" = tableau_colors[4], 
               "c5" = tableau_colors[5], 
               "t1" = tableau_colors[7], 
               "t3" = tableau_colors[8], 
               "Range" = tableau_colors[9])

uc$null <- as.character(uc$null)
uc[uc$null=="range", "null"] <- "Range"
uc$metric <- as.character(uc$metric)
uc[uc$metric=="fdis", "metric"] <- "FDis"
uc[uc$metric=="feve", "metric"] <- "FEve"
uc[uc$metric=="fric", "metric"] <- "FRic"
uc$null <- factor(uc$null, levels=c("compnet",
                                    "c1",
                                    "c2",
                                    "c3",
                                    "c4",
                                    "c5",
                                    "Range",
                                    "t1",
                                    "t3"))

ucplot <- ggplot(uc, aes(x=falsepos, y=truepos, color=null, fill=null, shape=metric))+
  xlim(c(0, 1))+
  ylim(c(0, 0.8))+
  coord_fixed()+
  geom_vline(xintercept=0.05, lty="dashed")+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=0)+
  geom_point(size=3)+
  scale_size_manual(values=c(6,2,2,2))+
  scale_shape_manual(values=c(15, 19, 17, 18), name="Metric", labels=c("compnet", "FDis", "FEve", "FRic"))+
  scale_color_manual(name="Model", values=nullmod_colors, labels=c("compnet",
                                                                   "C1",
                                                                   "C2",
                                                                   "C3",
                                                                   "C4",
                                                                   "C5",
                                                                   "Range",
                                                                   "T1",
                                                                   "T3"))+
  scale_fill_manual(name="Model", values=nullmod_colors, labels=c("compnet",
                                                                   "C1",
                                                                   "C2",
                                                                   "C3",
                                                                   "C4",
                                                                   "C5",
                                                                   "Range",
                                                                   "T1",
                                                                   "T3"))+
  xlab("False Detection Rate")+
  ylab("True Detection Rate")+
  theme_bw()
ucplot

corr$null <- as.character(corr$null)
corr[corr$null=="range", "null"] <- "Range"
corr$metric <- as.character(corr$metric)
corr[corr$metric=="fdis", "metric"] <- "FDis"
corr[corr$metric=="feve", "metric"] <- "FEve"
corr[corr$metric=="fric", "metric"] <- "FRic"
corr$null <- factor(corr$null, levels=c("compnet",
                                    "c1",
                                    "c2",
                                    "c3",
                                    "c4",
                                    "c5",
                                    "Range",
                                    "t1",
                                    "t3"))

corrplot <- ggplot(corr, aes(x=falsepos, y=truepos, color=null, fill=null, shape=metric))+
  xlim(c(0, 1))+
  ylim(c(0, 0.8))+
  coord_fixed()+
  geom_vline(xintercept=0.05, lty="dashed")+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=0)+
  geom_point(size=3)+
  scale_size_manual(values=c(6,2,2,2))+
  scale_shape_manual(values=c(15, 19, 17, 18), name="Metric", labels=c("compnet", "FDis", "FEve", "FRic"))+
  scale_color_manual(name="Model", values=nullmod_colors, labels=c("compnet",
                                                                   "C1",
                                                                   "C2",
                                                                   "C3",
                                                                   "C4",
                                                                   "C5",
                                                                   "Range",
                                                                   "T1",
                                                                   "T3"))+
  scale_fill_manual(name="Model", values=nullmod_colors, labels=c("compnet",
                                                                  "C1",
                                                                  "C2",
                                                                  "C3",
                                                                  "C4",
                                                                  "C5",
                                                                  "Range",
                                                                  "T1",
                                                                  "T3"))+
  xlab("False Detection Rate")+
  ylab("True Detection Rate")+
  theme_bw()+
  theme(axis.title.y=element_blank())
corrplot



ucpa$null <- as.character(ucpa$null)
ucpa[ucpa$null=="range", "null"] <- "Range"
ucpa$metric <- as.character(ucpa$metric)
ucpa[ucpa$metric=="fdis", "metric"] <- "FDis"
ucpa[ucpa$metric=="feve", "metric"] <- "FEve"
ucpa[ucpa$metric=="fric", "metric"] <- "FRic"
ucpa$null <- factor(ucpa$null, levels=c("compnet",
                                    "c1",
                                    "c2",
                                    "c3",
                                    "c4",
                                    "c5",
                                    "Range",
                                    "t1",
                                    "t3"))

ucpaplot <- ggplot(ucpa, aes(x=falsepos, y=truepos, color=null, fill=null, shape=metric))+
  xlim(c(0, 1))+
  ylim(c(0, 0.8))+
  coord_fixed()+
  ggtitle("Uncorrelated Traits")+
  geom_vline(xintercept=0.05, lty="dashed")+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=0)+
  geom_point(size=3)+
  scale_size_manual(values=c(6,2,2,2))+
  scale_shape_manual(values=c(15, 19, 17, 18), name="Metric", labels=c("compnet", "FDis", "FEve", "FRic"))+
  scale_color_manual(name="Model", values=nullmod_colors, labels=c("compnet",
                                                                   "C1",
                                                                   "C2",
                                                                   "C3",
                                                                   "C4",
                                                                   "C5",
                                                                   "Range",
                                                                   "T1",
                                                                   "T3"))+
  scale_fill_manual(name="Model", values=nullmod_colors, labels=c("compnet",
                                                                  "C1",
                                                                  "C2",
                                                                  "C3",
                                                                  "C4",
                                                                  "C5",
                                                                  "Range",
                                                                  "T1",
                                                                  "T3"))+
  xlab("False Detection Rate")+
  ylab("True Detection Rate")+
  theme_bw()+
  theme(axis.title.x=element_blank())
ucpaplot


corrpa$null <- as.character(corrpa$null)
corrpa[corrpa$null=="range", "null"] <- "Range"
corrpa$metric <- as.character(corrpa$metric)
corrpa[corrpa$metric=="fdis", "metric"] <- "FDis"
corrpa[corrpa$metric=="feve", "metric"] <- "FEve"
corrpa[corrpa$metric=="fric", "metric"] <- "FRic"
corrpa$null <- factor(corrpa$null, levels=c("compnet",
                                        "c1",
                                        "c2",
                                        "c3",
                                        "c4",
                                        "c5",
                                        "Range",
                                        "t1",
                                        "t3"))

corrpaplot <- ggplot(corrpa, aes(x=falsepos, y=truepos, color=null, fill=null, shape=metric))+
  xlim(c(0, 1))+
  ylim(c(0, 0.8))+
  coord_fixed()+
  ggtitle("Correlated Traits")+
  geom_vline(xintercept=0.05, lty="dashed")+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=0)+
  geom_point(size=3)+
  scale_size_manual(values=c(6,2,2,2))+
  scale_shape_manual(values=c(15, 19, 17, 18), name="Metric", labels=c("compnet", "FDis", "FEve", "FRic"))+
  scale_color_manual(name="Model", values=nullmod_colors, labels=c("compnet",
                                                                   "C1",
                                                                   "C2",
                                                                   "C3",
                                                                   "C4",
                                                                   "C5",
                                                                   "Range",
                                                                   "T1",
                                                                   "T3"))+
  scale_fill_manual(name="Model", values=nullmod_colors, labels=c("compnet",
                                                                  "C1",
                                                                  "C2",
                                                                  "C3",
                                                                  "C4",
                                                                  "C5",
                                                                  "Range",
                                                                  "T1",
                                                                  "T3"))+
  xlab("False Detection Rate")+
  ylab("True Detection Rate")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())
corrpaplot

summ_full <- (ucpaplot+corrpaplot)/(ucplot+corrplot)+plot_layout(guides = "collect")
summ_full
ggsave("summ_full.png", height=4.8, width=6)







ucplot <- ggplot(uc, aes(x=falsepos, y=truepos, color=null, fill=null, shape=metric))+
  ggtitle("Uncorrelated Traits")+
  xlim(c(0, 1))+
  ylim(c(0, 0.8))+
  coord_fixed()+
  geom_vline(xintercept=0.05, lty="dashed")+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=0)+
  geom_point(size=3)+
  scale_size_manual(values=c(6,2,2,2))+
  scale_shape_manual(values=c(15, 19, 17, 18), name="Metric", labels=c("compnet", "FDis", "FEve", "FRic"))+
  scale_color_manual(name="Model", values=nullmod_colors)+
  scale_fill_manual(name="Model", values=nullmod_colors)+
  xlab("False Detection Rate")+
  ylab("True Detection Rate")+
  theme_bw()
ucplot

corrplot <- ggplot(corr, aes(x=falsepos, y=truepos, color=null, fill=null, shape=metric))+
  ggtitle("Correlated Traits")+
  xlim(c(0, 1))+
  ylim(c(0, 0.8))+
  coord_fixed()+
  geom_vline(xintercept=0.05, lty="dashed")+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=0)+
  geom_point(size=3)+
  scale_size_manual(values=c(6,2,2,2))+
  scale_shape_manual(values=c(15, 19, 17, 18), name="Metric", labels=c("compnet", "FDis", "FEve", "FRic"))+
  scale_color_manual(name="Model", values=nullmod_colors)+
  scale_fill_manual(name="Model", values=nullmod_colors)+
  xlab("False Detection Rate")+
  ylab("True Detection Rate")+
  theme_bw()+
  theme(axis.title.y=element_blank())
corrplot

summ_abun <- ucplot+corrplot+plot_layout(guides = "collect")
summ_abun
ggsave("summ_abun.png", height=4, width=6)








summ_sub <- subset(summ, falsepos<0.1)
summ_sub_uc <- subset(summ_sub, version=="uc" | version=="ucpa")
summ_sub_corr <- subset(summ_sub, version=="corr" | version=="corrpa")

method_colors <- c("compnet" = tableau_colors[6], 
                    "fdis_c1" = tableau_colors[1], 
                    "fdis_c3" = tableau_colors[3], 
                    "fdis_c5" = tableau_colors[5])

method_shapes <- c("compnet" = 21, 
                   "fdis_c1" = 22, 
                   "fdis_c3" = 22, 
                   "fdis_c5" = 22)

method_names <- c("compnet"="compnet",
                  "fdis_c1"="FDis c1",
                  "fdis_c3"="FDis c3",
                  "fdis_c5"="FDis c5")

ucplot_sub <- ggplot(summ_sub_uc, aes(x=falsepos, y=truepos, color=methodname, fill=methodname, shape=methodname))+
  ylim(c(0, 0.5))+
  ggtitle("Uncorrelated Traits")+
  geom_vline(xintercept=0.05, lty="dashed")+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=0)+
  geom_point(size=6)+
  scale_shape_manual(name="Method", values=method_shapes, labels=method_names)+
  scale_color_manual(name="Method", values=method_colors, labels=method_names)+
  scale_fill_manual(name="Method", values=method_colors, labels=method_names)+
  xlab("False Detection Rate")+
  ylab("True Detection Rate")+
  theme_bw()+
  theme(aspect.ratio=1)+
  guides(size="none", alpha="none")
ucplot_sub



corrplot_sub <- ggplot(summ_sub_corr, aes(x=falsepos, y=truepos, color=methodname, fill=methodname, shape=methodname))+
  ylim(c(0, 0.5))+
  ggtitle("Correlated Traits")+
  geom_vline(xintercept=0.05, lty="dashed")+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=0)+
  geom_point(size=6)+
  scale_shape_manual(name="Method", values=method_shapes, labels=method_names)+
  scale_color_manual(name="Method", values=method_colors, labels=method_names)+
  scale_fill_manual(name="Method", values=method_colors, labels=method_names)+
  xlab("False Detection Rate")+
  ylab("True Detection Rate")+
  theme_bw()+
  theme(aspect.ratio=1)+
  guides(size="none", alpha="none")
corrplot_sub

ucplot_sub+corrplot_sub+plot_layout(guides = "collect")
ggsave("summ_sub.png", height=4, width=6)



