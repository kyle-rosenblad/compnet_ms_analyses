fd1site <- function(tvec, pavec){
  tvec <- sort(tvec[which(pavec==1),])
  EW <- diff(tvec)/2
  PEW <- EW/sum(EW)
  PEWalt <- PEW
  PEWalt[PEWalt>1/length(EW)] <- 1/length(EW)
  feve_tmp <- (sum(PEWalt)-(1/length(EW)))/(1-(1/length(EW)))
  
  fdis_tmp <- mean(abs(tvec-mean(tvec)))
  
  fric_tmp <- max(tvec)-min(tvec)
  
  cbind(fric_tmp, fdis_tmp, feve_tmp)
}

fdfull <- function(traitvector, pamat) {
  traitvector <- scale(traitvector)
  out <- t(apply(pamat, 1, function(row) fd1site(traitvector, row)))
  colnames(out) <- c("fric", "fdis", "feve")
  out
}





sum_cons_pairs <- function(x){
  if(length(x)<2) {
    return(numeric(0))
  }
  x[-1] + x[-length(x)]
}

fd1site_ab <- function(tvec, abvec){
  names(tvec) <- names(abvec)
  tvec <- tvec[which(abvec>0)]
  abvec <- abvec[which(abvec>0)]
  tvec <- sort(tvec)
  abvec <- abvec[names(tvec)]
  
  traitspaces <- diff(tvec)
  absums <- sum_cons_pairs(abvec)
  
  EW <- traitspaces/absums
  PEW <- EW/sum(EW)
  PEWalt <- PEW
  PEWalt[PEWalt>1/length(EW)] <- 1/length(EW)
  feve_tmp <- (sum(PEWalt)-(1/length(EW)))/(1-(1/length(EW)))
  
  traitcentroid <- sum(abvec*tvec)/sum(abvec)
  centroiddists <- abs(tvec-traitcentroid)
  
  fdis_tmp <- sum(abvec*centroiddists)/sum(abvec)
  
  fric_tmp <- max(tvec)-min(tvec)
  
  cbind(fric_tmp, fdis_tmp, feve_tmp)
}

fdfull_ab <- function(traitvector, abmat) {
  traitvector <- scale(traitvector)
  out <- t(apply(abmat, 1, function(row) fd1site_ab(traitvector, row)))
  colnames(out) <- c("fric", "fdis", "feve")
  out
}
