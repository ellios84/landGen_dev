# Function to get Cross-Entropy decay based on sNMF analysis 

plotCeDecay <- function(snmf_res = NULL, maxK = NULL, figLab = NULL) {
  
  # Creating a "cross-entropy" object for each K 
  for (i in 1:maxK) assign(paste0("ce_k", sprintf("%02d", i)), LEA::cross.entropy(snmf_res, K = i))
  
  # Individuating these objects in the working space
  tmp <- ls()[grep(ls(), pattern = "ce_k")]
  
  # initialising the medians and notches arrays
  medians <- array(data = NA, dim = maxK)
  notches <- array(data = NA, dim = maxK)
  
  for (i in 1:maxK) {
    
    # Updating the ce_ki objects with correct row names (one run per line) and a column specifying the right K
    ce_ki <- as.data.frame(get(tmp[i]))
    rownames(ce_ki) <- paste0("K", sprintf("%02d", i), "_run", sprintf("%03d", 1:nrow(ce_ki)))
    ki <- colnames(ce_ki)
    ki <- substr(ki, 5 , nchar(ki))
    ce_ki$V2 <- ki
    colnames(ce_ki) <- c("CE", "K")
    assign(paste0("ce_k", sprintf("%02d", i)), ce_ki)
    
    # Calculating the notch of the distribution
    medians[i] <- median(ce_ki$CE)
    notches[i] <- 1.58 * (IQR(ce_ki$CE)/sqrt(nrow(ce_ki)))
    
    rm(ce_ki, ki)
  }
  
  # Create a unique data frame containing all results
  ce <- data.frame()
  for (i in 1:maxK) {
    ce <- rbind.data.frame(ce, get(tmp[i]))
  }
  rm(i, tmp)
  ce$K <- as.numeric(ce$K)
  
  # Finding overlaps (if any) between notches
  medCE <- by(data = ce$CE, INDICES = ce$K, FUN = median)
  Klab <- paste0("K", names(medCE))
  medCE <- as.vector(medCE)
  names(medCE) <- Klab
  bestK <- which.min(medCE)
  notch_df <- data.frame(L = medians - notches, U = medians + notches)
  best_notch <- data.frame(L = medians[bestK] - notches[bestK], U = medians[bestK] + notches[bestK])
  notch_results <- array(data = NA, dim = nrow(notch_df))
  for (i in 1:length(notch_results)) {
    if(notch_df[i, 1] <= best_notch[1, 2]) {# & notch_df[i, 1] >= best_notch[1, 1]) {
      notch_results[i] <- 1 
    } else notch_results[i] <- 0 
  }; rm(i)
  names(notch_results) <- Klab
  notch_results <- which(notch_results == 1)
  
  # Plot Ce decay results
  jpeg(paste0(figLab, "_Ce_decay.jpeg", sep = ""), width = 8, height = 6, units = 'in', res = 800)
  par(mar = c(4, 4, 1, 1), oma=c(2, 1, 2, 1))
  plot(range(ce$K), range(ce$CE), t = "n", xlab = "Number of ancestral populations", ylab = "Cross-Entropy", ylim = range(ce$CE), xlim = range(ce$K), las = 1, cex.axis = 0.8, cex.lab = 1.3)
  points(CE~jitter(x = K, factor = 1), data = ce, pch = 16, cex = .3, col = "gray")
  bxplt_col <- rep(scales::alpha("gray92", .5), maxK)
  border_col <- rep("gray92", maxK)
  if (length(notch_results) > 0) {
    bxplt_col[notch_results] <- scales::alpha("red", .5)
    border_col[notch_results] <- "red" 
  }
  boxplot(CE~K, data = ce, add = TRUE, xaxt = "n", yaxt = "n", border = border_col,
          col = bxplt_col, pch = 16, cex = 1, outline = FALSE, notch = TRUE)
  
  # Highlighting of the proposed K(s)
  points(bestK, medCE[bestK], pch = 21, bg = "red", cex = .7)
  dev.off()
  
  # Suggested K
  cat(paste0("Minimum median cross-entropy (MMCE) at K = ", bestK, " (MMCE = ", round(medCE[bestK], 4), ")"))
  cat("\n")
  
  # Minimum Cross-Entropy (MCE)
  MCE <- ce[which.min(ce$CE), ]
  cat(paste0("Minimum absolute cross-entropy (MACE) at K = ", MCE[2], ", run = ", substr(x = rownames(MCE), start = 8, stop = nchar(rownames(MCE))), " (MACE = ", round(MCE[1], 4), ")"))   
  cat("\n")
  
  # returns Ce decay results
  list.res <- list()
  list.res[[1]] <- ce
  list.res[[2]] <- notch_df
  names(list.res) <- c("Cross-Entropy", "CI95%")
  return(list.res)
  
}
