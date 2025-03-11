options(scipen=999)
library(foreach)
library(doParallel)
registerDoParallel(cores=8)

# Set up Thompson et al.'s (2020) parameter space
param_levels <- readRDS("param_levels.RDS")

# Set the number of species
species <- readRDS("species.RDS")

sim_com_output <- readRDS("sim_com_output.RDS")

set.seed(2347)
compnet_res <- foreach(g=1:nrow(param_levels)) %dopar% {
  gc()
  library(compnet)
  mat <- sim_com_output[[g]][[1]]
  mat[mat>0] <- 1
  sums <- rowSums(mat)
  coocc_sites <- length(subset(sums, sums>1))
  rich <- colSums(mat)
  rich[rich>1] <- 1
  rich <- sum(rich)
  
  if(coocc_sites==0 | rich<3){
    result <- "insufficient"
  }
  
  if(coocc_sites>0 & rich>2){
    sp_df <- sim_com_output[[g]][[2]]
    colnames(mat) <- sp_df$sp
    mat <- mat[, colSums(mat)>0]
    spp <- colnames(mat)
    sp_df <- sp_df[sp_df$sp%in%spp, ]
    rownames(sp_df) <- sp_df$sp
    
    mod <- buildcompnet(presabs = mat,
                        spvars_multi_int = sp_df[c("ndtrait", "domtrait")],
                        rank=1,
			prior_betas_scale=1)
    result <- mod$stanmod_samp$beta_dy[,1]
  }
  result
}
saveRDS(compnet_res, "compnet_results.RDS")
