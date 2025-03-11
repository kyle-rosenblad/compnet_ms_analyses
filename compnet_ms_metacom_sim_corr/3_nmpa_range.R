options(scipen=999)
library(foreach)
library(doParallel)
registerDoParallel(cores=20) # adjust number of cores for different environment

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
  
  envnum <- param_levels[g, "rep"]
  randenv <- readRDS(paste("thompson_env_", envnum, ".RDS", sep=""))
  randenv <- subset(randenv, time==max(randenv$time))
  randenv$site <- paste(randenv$x, randenv$y, "1", sep="_")
  randenv <- subset(randenv, site%in%rownames(mat))
  env_data <- c(randenv$envx_t)
  names(env_data) <- randenv$site

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
    
    matsumms <- data.frame(site=rownames(mat), rich=rowSums(mat))
    sp_df2 <- sp_df[c("ndtrait")]
    rownames(sp_df2) <- sp_df$sp
    
    source("fd_functions.R")
    obsfd <- fdfull_ab(sp_df2, mat)
    obsfd <- as.data.frame(obsfd)
    names(obsfd) <- c("fric", "fdis", "feve")
    matsumms <- cbind(matsumms, obsfd)
    
    # Function to randomize the abundance matrix
    randomize_matrix <- function(ab_matrix) {
      # Number of sites and species
      n_sites <- nrow(ab_matrix)
      n_species <- ncol(ab_matrix)
      
      # Create a new matrix to store randomized data
      randomized_matrix <- ab_matrix
      
      # Iterate over each site
      for (site in 1:n_sites) {
        # Get the environmental value for the current site
        site_env_value <- env_data[site]
        
        # Determine which species are in range of this site's environmental value
        species_in_range <- which(
          apply(ab_matrix, 2, function(species_col) {
            species_sites <- which(species_col == 1)
            range_min <- min(env_data[species_sites])
            range_max <- max(env_data[species_sites])
            range_min <= site_env_value && site_env_value <= range_max
          })
        )
        
        # Extract the subset of the presence-absence data for these species
        ab_subset <- ab_matrix[site, species_in_range]
        
        # Resample (reshuffle) within this subset
        if(length(ab_subset)>1){
          randomized_subset <- sample(ab_subset)
        }
        if(length(ab_subset)<2){
          randomized_subset <- ab_subset
        }
        
        # Place the randomized subset back into the matrix
        randomized_matrix[site, species_in_range] <- randomized_subset
      }
      
      colnames(randomized_matrix) <- colnames(mat)
      # Remove species that are absent from all sites after randomization
      col_sums <- colSums(randomized_matrix)
      randomized_matrix <- randomized_matrix[, col_sums > 0, drop = FALSE]
      
      return(randomized_matrix)
    }
    randmats <- lapply(1:999, function(i) randomize_matrix(mat))
    
    compute_fd_random <- function(inmat){
      sp_df3 <- data.frame(ndtrait=sp_df2[rownames(sp_df2)%in%colnames(inmat),])
      rownames(sp_df3) <- sp_df[sp_df$sp%in%colnames(inmat), "sp"]
      fdfull_ab(sp_df3, inmat)
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
saveRDS(nmpa_res, "nmpa_results_range.RDS")
