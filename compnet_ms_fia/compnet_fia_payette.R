options(scipen=999)
library(sf)
library(rtry)

# get species for which try data were requested
try_species <- read.csv("try_species.csv")

# read and process western US tree data from Rosenblad et al. (2023)
sites <- readRDS("css2.RDS")
trees <- readRDS("data2.RDS")
trees <- subset(trees, LAT_LON_SUBP%in%sites$LAT_LON_SUBP)
trees$binomial <- paste(trees$GENUS, trees$SPECIES, sep=" ")

# get species dataframe for fia data
species <- unique(trees[c("GENUS", "SPECIES", "SPCD")])
species$binomial <- paste(species$GENUS, species$SPECIES, sep=" ")
species <- merge(species, try_species, by.x="binomial", by.y="AccSpeciesName",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)

# enter TRY IDs for species that didn't get one in the merge
subset(species, is.na(AccSpeciesID))
species <- subset(species, binomial!="Tree evergreen")
subset(species, is.na(AccSpeciesID))
species[species$binomial=="Pinus washoensis", "AccSpeciesID"] <- 446578
species[species$binomial=="Cupressus arizonica", "AccSpeciesID"] <- 429665
species[species$binomial=="Abies shastensis", "AccSpeciesID"] <- 200049
species[species$binomial=="Lithocarpus densiflorus", "AccSpeciesID"] <- 244314
species[species$binomial=="Cupressus macrocarpa", "AccSpeciesID"] <- 429674
species[species$binomial=="Chamaecyparis nootkatensis", "AccSpeciesID"] <- 412055
species[species$binomial=="Cupressus sargentii", "AccSpeciesID"] <- 429677

# read try data
trydata <- rtry_import("trydata.txt")
trydata <- rtry_exclude(trydata, (DataID %in% 413) & (OrigValueStr %in% c("juvenile", "saplings")), baseOn = ObservationID)

# process height data, convert to common units
ht <- subset(trydata, TraitID==3106)
table(ht$OrigUnitStr) # convert all to m
ht$val <- as.numeric(ht$OrigValueStr)
ht[ht$OrigUnitStr=="cm", "val"] <- ht[ht$OrigUnitStr=="cm", "val"]/100
ht[ht$OrigUnitStr=="feet", "val"] <- ht[ht$OrigUnitStr=="feet", "val"]*0.3048
ht <- tapply(ht$val, INDEX=ht$AccSpeciesID, FUN=mean, na.rm=TRUE)
ht <- data.frame(AccSpeciesID=names(ht), ht=ht)

# process wood density data, convert to common units
wd <- subset(trydata, TraitID==4)
table(wd$OrigUnitStr) # convert all to g / cm^3
wd$val <- as.numeric(wd$OrigValueStr)
wd[wd$OrigUnitStr=="g/dm3", "val"] <- 0.001*wd[wd$OrigUnitStr=="g/dm3", "val"]
wd[wd$OrigUnitStr=="kg/m3", "val"] <- 0.001*wd[wd$OrigUnitStr=="kg/m3", "val"]
wd[wd$OrigUnitStr=="mg/cm3", "val"] <- 0.001*wd[wd$OrigUnitStr=="mg/cm3", "val"]
wd <- subset(wd, OrigUnitStr!="t/m3")
wd <- tapply(wd$val, INDEX=wd$AccSpeciesID, FUN=mean, na.rm=TRUE)
wd <- data.frame(AccSpeciesID=names(wd), wd=wd)

# combine traits into one data set
trydata2 <- merge(ht, wd, by="AccSpeciesID", all.x=TRUE, all.y=TRUE)

# filter down to only data from payette national forest
nf <- st_read(dsn = getwd(), layer = "S_USA.AdministrativeForest")
trees <- st_as_sf(trees, coords=c("LON", "LAT"), crs=st_crs(nf))
sort(nf$FORESTNAME)
sf_use_s2(FALSE)
nf2 <- subset(nf, FORESTNAME=="Payette National Forest")
trees2 <- st_filter(trees, nf2)
sites <- subset(sites, LAT_LON%in%trees2$LAT_LON)
trees2$binomial <- paste(trees2$GENUS, trees2$SPECIES, sep=" ")
species <- subset(species, SPCD%in%trees2$SPCD)
species <- merge(species, trydata2, by="AccSpeciesID",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)
trees2 <- merge(trees2, species[c("binomial", "AccSpeciesID")], by="binomial",
                all.x=TRUE, all.y=FALSE, sort=FALSE)
presabs <- matrix(NA, nrow=length(unique(sites$LAT_LON_SUBP)), ncol=nrow(species))
rownames(presabs) <- unique(sites$LAT_LON_SUBP)
colnames(presabs) <- species$SPCD
trees2 <- subset(trees2, INVYR==survey_2_year)

# set up presence-absence matrix for compnet analysis
for(i in 1:ncol(presabs)){
  sites_tmp <- unique(subset(trees2, SPCD==colnames(presabs)[i])$LAT_LON_SUBP)
  presabs[,i] <- as.numeric(rownames(presabs)%in%sites_tmp)
}

# set up trait matrix for compnet analysis
traitdata <- species[c("SPCD", "binomial", "ht", "wd")]
rownames(traitdata) <- traitdata$SPCD
traitdata$loght <- log10(traitdata$ht)
traitdata
cor(traitdata[c("loght", "wd")])

## compnet analyses
library(compnet)
library(DHARMa)
set.seed(8316)
mod <- buildcompnet(presabs=presabs,
                    spvars_multi_int=traitdata[c("loght", "wd")],
                    warmup=1000,
                    iter=2000,
                    rank=1,
                    prior_betas_scale=1,
                    adapt_delta=0.99)
mod2 <- buildcompnet(presabs=presabs,
                    spvars_multi_int=traitdata[c("loght")],
                    warmup=1000,
                    iter=2000,
                    rank=1,
                    prior_betas_scale=1,
                    adapt_delta=0.99)
summarize_compnet(mod)
summarize_compnet(mod2)

ppred <- postpredsamp(mod)
fpr <- apply(ppred, 1, mean)
modcheck <- createDHARMa(simulatedResponse = ppred,
                         observedResponse = mod$d$both,
                         fittedPredictedResponse = fpr,
                         integerResponse = TRUE)
testDispersion(modcheck)
testQuantiles(modcheck)
testUniformity(modcheck)
testZeroInflation(modcheck)
gofstats(mod)




histsdf1 <- data.frame(es=mod$stanmod_samp$beta_dy[,1],
                       model="Full")
histsdf2 <- data.frame(es=mod2$stanmod_samp$beta_dy[,1],
                       model="Single-Trait")
histsdf <- rbind(histsdf1, histsdf2)

sp <- scatter_interaction(mod, "loght",
                          xlabel="Log(Height [m])",
                          ci_width=0.95, thin=FALSE,
                          ymax=1)
sp_wd <- scatter_interaction(mod, "wd",
                              xlabel="Wood Density [g/mL]",
                          ci_width=0.95, thin=FALSE,
                          ymax=1)
sp2 <- scatter_interaction(mod2, "loght",
                           xlabel="Log(Height [m])",
                           ci_width=0.95, thin=FALSE,
                           ymax=1)

library(ggplot2)
sp <- sp+
  ggtitle("Full Model")
sp_wd <- sp_wd+
  ggtitle("Full Model")
sp2 <- sp2+
  ggtitle("Single-Trait Model")


eshists <- ggplot(histsdf, aes(x=es, group=model, color=model))+
  scale_color_manual(values=c("red", "blue"),
                     name="Model")+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_density(linewidth=1)+
  xlab("Standardized Height Interaction Effect")+
  ylab("Density")+
  theme_bw()+
  theme(aspect.ratio=1)
eshists

library(patchwork)
sp+sp_wd+sp2+eshists
ggsave("payette.png", width=9, height=6)
