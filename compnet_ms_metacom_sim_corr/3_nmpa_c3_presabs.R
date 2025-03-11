options(scipen=999)
library(foreach)
library(doParallel)
#registerDoParallel(cores=20)
registerDoParallel(cores=8)

# Set up Thompson et al.'s (2020) parameter space
param_levels <- readRDS("param_levels.RDS")

# Set the number of species
species <- readRDS("species.RDS")

sim_com_output <- readRDS("sim_com_output.RDS")

set.seed(275)
nmpa_res <- foreach(g=1:nrow(param_levels), .combine="rbind") %dopar% {
  gc()
  library(vegan)
  mat <- sim_com_output[[g]][[1]]
  mat2 <- mat
  mat2[mat2>0] <- 1
  sums <- rowSums(mat2)
  coocc_sites <- length(subset(sums, sums>1))
  rich <- colSums(mat2)
  rich[rich>1] <- 1
  rich <- sum(rich)
  
  if(coocc_sites==0 | rich<3){
    p_final_fric <- "insufficient"
    p_final_fdis <- "insufficient"
    p_final_feve <- "insufficient"
  }
  
  if(coocc_sites>0 & rich>=3){
    sp_df <- sim_com_output[[g]][[2]]
    colnames(mat) <- sp_df$sp
    mat <- mat[, colSums(mat)>0]
    spp <- colnames(mat)
    sp_df <- sp_df[sp_df$sp%in%spp, ]
    sp_df$spid <- 1:nrow(sp_df)
    
    mat <- mat[rowSums(mat)>0,]
    
    mat[mat>0] <- 1
    
    matsumms <- data.frame(site=rownames(mat), rich=rowSums(mat))
    sp_df2 <- sp_df[c("ndtrait")]
    rownames(sp_df2) <- sp_df$sp
    
    source("fd_functions.R")
    obsfd <- fdfull(sp_df2, mat)
    obsfd <- as.data.frame(obsfd)
    names(obsfd) <- c("fric", "fdis", "feve")
    matsumms <- cbind(matsumms, obsfd)
    
    reshuffle_nonzero <- function(vec) {
      nonzero_indices <- which(vec != 0)
      nonzero_values <- vec[nonzero_indices]
      if(length(nonzero_values)>1){
        shuffled_values <- sample(nonzero_values)
        reshuffled_vec <- vec
        reshuffled_vec[nonzero_indices] <- shuffled_values
      }
      if(length(nonzero_values)<2){
        reshuffled_vec <- vec
      }
      return(reshuffled_vec)
    }
    reshuffle_cols <- function(inmat) {
      apply(inmat, 2, reshuffle_nonzero)
    }
    fix_col_names <- function(inmat){
      colnames(inmat) <- colnames(mat)
      inmat
    }
    randmats <- lapply(1:999, function(x) reshuffle_cols(mat))
    randmats <- lapply(randmats, fix_col_names)
    
    compute_fd_random <- function(inmat){
      fdfull(sp_df2, inmat)
    }
    
    process_random_matrices <- function(mat_list) {
      result_list <- lapply(seq_along(mat_list), function(idx) {
        output_matrix <- compute_fd_random(mat_list[[idx]])
        output_matrix <- cbind(output_matrix, rep = idx)
        return(output_matrix)
      })
      
      combined_result <- do.call(rbind, result_list)
      return(combined_result)
    }
    
    randmatsumms <- process_random_matrices(randmats)
    sitesvec <- rownames(randmatsumms)
    randmatsumms <- as.data.frame(randmatsumms)
    names(randmatsumms) <- c("fric", "fdis", "feve", "rep")
    randmatsumms$site <- sitesvec
    
    for(l in 1:nrow(matsumms)){
      if(!is.na(matsumms[l,"fdis"])){
        summstmp <- subset(randmatsumms, site==matsumms[l,"site"])
        matsumms[l, "fric_z"] <- (matsumms[l, "fric"]-mean(summstmp$fric))/sd(summstmp$fric)
        matsumms[l, "fdis_z"] <- (matsumms[l, "fdis"]-mean(summstmp$fdis))/sd(summstmp$fdis)
        matsumms[l, "feve_z"] <- (matsumms[l, "feve"]-mean(summstmp$feve))/sd(summstmp$feve)
      }
    }
    
    tmp <- subset(matsumms, !is.na(fric_z))
    if(nrow(tmp)>0){
      wtest <- wilcox.test(tmp$fric_z, alternative="greater")
      p_final_fric <- wtest$p.value
      fric_z_bar <- mean(tmp$fric_z)
    }
    if(nrow(tmp)==0){
      p_final_fric <- "insufficient"
    }
    
    tmp <- subset(matsumms, !is.na(fdis_z))
    if(nrow(tmp)>0){
      wtest <- wilcox.test(tmp$fdis_z, alternative="greater")
      p_final_fdis <- wtest$p.value
      fdis_z_bar <- mean(tmp$fdis_z)
    }
    if(nrow(tmp)==0){
      p_final_fdis <- "insufficient"
    }
    
    tmp <- subset(matsumms, !is.na(feve_z))
    if(nrow(tmp)>0){
      wtest <- wilcox.test(tmp$feve_z, alternative="greater")
      p_final_feve <- wtest$p.value
      feve_z_bar <- mean(tmp$feve_z)
    }
    if(nrow(tmp)==0){
      p_final_feve <- "insufficient"
    }
  }
  c(p_final_fric, p_final_fdis, p_final_feve)
}
nmpa_res <- as.data.frame(nmpa_res)
names(nmpa_res) <- c("p_fric", "p_fdis", "p_feve")
nmpa_res$p_fric <- as.numeric(nmpa_res$p_fric)
nmpa_res$p_fdis <- as.numeric(nmpa_res$p_fdis)
nmpa_res$p_feve <- as.numeric(nmpa_res$p_feve)
saveRDS(nmpa_res, "nmpa_results_c3_presabs.RDS")
