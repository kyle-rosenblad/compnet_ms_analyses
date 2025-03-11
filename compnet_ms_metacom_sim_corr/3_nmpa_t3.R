options(scipen=999)
library(foreach)
library(doParallel)
registerDoParallel(cores=20)

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
    
    matsumms <- data.frame(site=rownames(mat), rich=rowSums(mat))
    sp_df2 <- sp_df[c("ndtrait")]
    rownames(sp_df2) <- sp_df$sp
    
    source("fd_functions.R")
    obsfd <- fdfull_ab(sp_df2, mat)
    obsfd <- as.data.frame(obsfd)
    names(obsfd) <- c("fric", "fdis", "feve")
    matsumms <- cbind(matsumms, obsfd)
    
    spabun <- colSums(mat)
    spabun_rare <- spabun[spabun<=median(spabun)]
    spabun_common <- spabun[spabun>median(spabun)]
    raresp <- subset(sp_df, sp%in%names(spabun_rare))
    commonsp <- subset(sp_df, sp%in%names(spabun_common))
    
    reshuffle_traits <- function() {
      raresp_tmp <- raresp
      commonsp_tmp <- commonsp
      if(nrow(raresp_tmp)>1){
        raresp_tmp$ndtrait <- sample(raresp_tmp$ndtrait)
      }
      if(nrow(raresp_tmp)<2){
        raresp_tmp$ndtrait <- raresp_tmp$ndtrait
      }
      if(nrow(commonsp_tmp)>1){
        commonsp_tmp$ndtrait <- sample(commonsp_tmp$ndtrait)
      }
      if(nrow(commonsp_tmp)<2){
        commonsp_tmp$ndtrait <- commonsp_tmp$ndtrait
      }
      outdf <- rbind(raresp_tmp, commonsp_tmp)
      outdf <- outdf[match(sp_df$sp, outdf$sp),]
      rownames(outdf) <- outdf$sp
      return(outdf[c("ndtrait")])
    }
    
    randtraits <- lapply(1:999, function(x) reshuffle_traits())

    compute_fd_random <- function(indf){
      fdfull_ab(indf, mat)
    }
    
    process_random_traits <- function(df_list) {
      result_list <- lapply(seq_along(df_list), function(idx) {
        output_matrix <- compute_fd_random(df_list[[idx]])
        output_matrix <- cbind(output_matrix, rep = idx)
        return(output_matrix)
      })
      
      combined_result <- do.call(rbind, result_list)
      return(combined_result)
    }
    
    randtraitsumms <- process_random_traits(randtraits)
    sitesvec <- rownames(randtraitsumms)
    randtraitsumms <- as.data.frame(randtraitsumms)
    names(randtraitsumms) <- c("fric", "fdis", "feve", "rep")
    randtraitsumms$site <- sitesvec
    
    for(l in 1:nrow(matsumms)){
      if(!is.na(matsumms[l,"fdis"])){
        summstmp <- subset(randtraitsumms, site==matsumms[l,"site"])
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
saveRDS(nmpa_res, "nmpa_results_t3.RDS")
