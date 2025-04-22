# idealmente lfmm_analysis --> gea(lfmm, rda, pRDA, varpart neutral vs adaptive)

lfmm_analysis <- function (Q = NULL, dist = "Aitchison", lfmm_in = NULL, E = NULL, 
                           minK = 1, maxK = NULL, analysisLab = "Ejemplo", figure = TRUE) {
  
  # Loading the Q-matrix
  Q <- read.table(file = Q, header = FALSE)
  
  if (dist == "Euclidean") {
    # Euclidean distance
    Qdist <- dist(Q, method = "euclidean") 
  }
  
  if (dist == "Aitchison") {
    # Q-dist: Aitchison distance (accounts for the relative scale property of compositional data)
    Qdist <- robCompositions::aDist(x = Q)
  }
  
  if (dist == "Bray") {
    # Bray-Curtis distance
    Qdist <- vegan::vegdist(x = Q, method = "bray")
  }
      
  # Loading the environmental matrix 
  E <- read.table(file = E, header = TRUE)
  
  # Loading genetic matrix for LFMMs
  lfmm_in <- data.table::fread(lfmm_in)
  
  # List where to store LFMM results
  lfmm_ls <- list()
  
  # Array to store for Mantel statistics (r)
  mantel_R <- array()
  
  # Array to store RDA adjusted R squared
  rda_R <- array()
  
  # Array to store pRDA adjusted R squared
  pRda_R <- array() 
  
  # my_qvalue function
  my_qvalue <- function(x) {
    qval <- qvalue::qvalue(x)
    qval <- qval$qvalues
    return(qval)
  }
  
  # Running the analysis... 
  for (i in minK:maxK) {
    
    # Lfmm
    lfmm_i <- lfmm::lfmm_ridge(Y = lfmm_in, X = E, K = i)
    pval_i <- lfmm::lfmm_test(Y = lfmm_in, X = E, lfmm = lfmm_i, calibrate = "gif")
    pval_i <- as.data.frame(pval_i$calibrated.pvalue)
    
    # q-value calculation
    pval_i <- apply(pval_i, 2, my_qvalue)
    pval_i <- as.data.frame(pval_i)
    lfmm_ls[[i]] <- pval_i
    
    # Latent factor scores
    ps_i <- lfmm_i$U

    # Mantel test
    LFdist_i <- dist(ps_i, method = "euclidean") # Euclidean distance
    mantel_i <- vegan::mantel(Qdist, LFdist_i, method = "spearman", permutations = 99, na.rm = TRUE, parallel = 4)
    mantel_R[i] <- mantel_i$statistic # pval[i] <- mantel_i$signif
    
    # RDA
    rda_i <- vegan::rda(Q ~ ps_i)
    rda_R[i] <- vegan::RsquareAdj(rda_i)$adj.r.squared

    # pRDA
    pRda_i <- vegan::rda(lfmm_in ~ ps_i + Condition(as.matrix(E)))
    pRda_R[i] <- vegan::RsquareAdj(pRda_i)$adj.r.squared
    
    cat(paste0("K=", i, " - Mantel statistic R: ", round(mantel_R[i], 4), "; Q var. explained by K: ", round(rda_R[i], 4), "; Molecular var. explained by K (conditioned on environment): ", round(pRda_R[i], 4)))
    cat("\n")
    
  }
  
  # Plot
  if (figure == TRUE) {
    jpeg(paste0(analysisLab, ".jpeg"), width = 8, height = 4.5, units = 'in', res = 800)
    par(mfrow = c(1, 2))
    plot(minK:maxK, mantel_R, las = 1, xlab = "Number of latent factors", ylab = "Mantel statistic r", pch = 16, ylim = c(-1, 1)); abline(h = 0, lty = 5)
    plot(minK:maxK, rda_R, las = 1, xlab = "Number of latent factors", ylab = "Adjusted R squared", pch = 16, ylim = c(0, 1))
    dev.off()
  }
  
  # Return (da capire come ritornare i q-val per ogni K)
  save(lfmm_ls, file = paste0(analysisLab, "_lfmmResults.RData"))
  res <- list(mantel_R, rda_R, pRda_R)
  names(res) <- c("mantel_R", "rda_R", "pRda_R")
  cat("\n")
  return(res)

}
