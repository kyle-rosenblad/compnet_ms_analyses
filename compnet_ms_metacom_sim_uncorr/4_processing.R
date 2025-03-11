options(scipen=999)
library(ggplot2)
library(ggthemes)
library(patchwork)

param_levels <- readRDS("param_levels.RDS")

r <- readRDS("compnet_results.RDS")
for(i in 1:length(r)){
  if(class(r[[i]])!="character"){
    param_levels[i, "compnet"] <- mean(r[[i]]<0)
  }
  if(class(r[[i]])=="character"){
    param_levels[i, "compnet"] <- NA
  }
}

r <- readRDS("compnet_results_neg.RDS")
for(i in 1:length(r)){
  if(class(r[[i]])!="character"){
    param_levels[i, "compnet_neg"] <- mean(r[[i]]<0)
  }
  if(class(r[[i]])=="character"){
    param_levels[i, "compnet_neg"] <- NA
  }
}
head(param_levels)
param_levels$index <- paste(param_levels$disp_val, param_levels$nb, sep="_")

nmpa_metrics <- c("fric", "fdis", "feve")
nmpa_nulls <- c("c1", "c2", "c3", "c4", "c5", "range", "t1", "t3")

for(i in 1:length(nmpa_nulls)){
  r <- readRDS(paste("nmpa_results_", nmpa_nulls[i], ".RDS", sep=""))
  r <- as.data.frame(r)
  names(r) <- paste(names(r), nmpa_nulls[i], sep="_")
  param_levels <- cbind(param_levels, r)
}
for(i in 1:length(nmpa_nulls)){
  r <- readRDS(paste("nmpa_results_", nmpa_nulls[i], "_neg.RDS", sep=""))
  r <- as.data.frame(r)
  names(r) <- paste(names(r), nmpa_nulls[i], "neg", sep="_")
  param_levels <- cbind(param_levels, r)
}
head(param_levels)

# are there any data sets for which compnet's result is NA, and another method's result is not?
table(c(as.matrix(subset(param_levels, is.na(compnet))[,7:30])), useNA="always")
# no

# what about for negative cases?
table(c(as.matrix(subset(param_levels, is.na(compnet_neg))[,31:54])), useNA="always")
# no

param_levels_pos <- subset(param_levels, !is.na(compnet))
param_levels_neg <- subset(param_levels, !is.na(compnet_neg))

nmpa_methods <- merge(nmpa_metrics, nmpa_nulls)
names(nmpa_methods) <- c("metric", "null")
nmpa_methods$methodname <- paste(nmpa_methods$metric, nmpa_methods$null, sep="_")

for(i in 1:nrow(nmpa_methods)){
  pval <- paste("p", nmpa_methods[i, "methodname"], sep="_")
  restmp <- c(as.matrix(param_levels_pos[pval]))
  nmpa_methods[i, "truepos"] <- sum(restmp<0.05 & !is.na(restmp))/length(restmp)
}
for(i in 1:nrow(nmpa_methods)){
  pval <- paste("p", nmpa_methods[i, "methodname"], "neg", sep="_")
  restmp <- c(as.matrix(param_levels_neg[pval]))
  nmpa_methods[i, "falsepos"] <- sum(restmp<0.05 | is.na(restmp))/length(restmp)
}
nmpa_methods

methods <- nmpa_methods
methods[nrow(methods)+1, "methodname"] <- "compnet"
methods[nrow(methods), "truepos"] <- mean(param_levels_pos$compnet>0.95)
methods[nrow(methods), "falsepos"] <- mean(param_levels_neg$compnet_neg>0.95)
methods[nrow(methods), "metric"] <- "compnet"
methods[nrow(methods), "null"] <- "compnet"
methods$null <- factor(methods$null, levels=c("compnet", "c1", "c2", "c3", "c4", "c5", "t1", "t3", "range"))

ggplot(methods, aes(x=falsepos, y=truepos, color=null, fill=null, shape=metric))+
  geom_point(size=3)+
  geom_vline(xintercept=0.05, lty="dashed")+
  scale_shape_manual(values=c(21:24), name="Metric")+
  scale_fill_tableau(name="Null Model")+
  scale_color_tableau(name="Null Model")+
  xlab("False Detection Rate")+
  ylab("True Detection Rate")+
  theme_bw()+
  theme(aspect.ratio=1)

saveRDS(methods, "../compnet_ms_metacom_sim_summary/summary_uncorr.RDS")
saveRDS(param_levels, "../compnet_ms_metacom_sim_summary/pl_uncorr.RDS")
