# IMPRTANT NOTE: in its current version, this funciton works with ped/map file extensions

pairwiseFst <- function(pedmap = NULL, path2plink2 = "./plink2.exe", coordInd = NULL) {
  
  # creating the '--within' file for pirwise Fst calculation (FIDs IDs FIDs)
  within <- data.table::fread(paste0(pedmap, ".ped"), header = FALSE)[, 1:2]
  within <- data.frame(FID = within$V1, ID = within$V2, FID = within$V1)
  write.table(x = within, file = "within.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)  
  
  # Calculate pairwise Fst (Plink 2.0)  
  system(paste0(path2plink2, " --pedmap ", pedmap, " --allow-extra-chr --within within.txt --fst CATPHENO method=hudson"))
  
  # Pairwise Fst distribution 
  pwFst <- read.table("plink2.fst.summary", header = FALSE)
  
  jpeg("pairwiseFst-boxplot.jpeg", width = 7, height = 7, units = 'in', res = 800)
  par(mar = c(1, 6, 1, 1), oma=c(2, 1, 2, 1))
  plot(rep(1, nrow(pwFst)), pwFst$V3, axes = FALSE, t = "n", las = 1, ylab = expression(paste("Pairwise ", italic(F)[ST])))
  points(jitter(x = rep(1, nrow(pwFst)), factor = 1.5, amount = .07), pwFst$V3, pch = 16, col = scales::alpha("lightgray", .5), yaxt = "n")
  boxplot(pwFst$V3, add = TRUE, border = "black", col = scales::alpha("lightgray", .5), outline = FALSE, las = 1)
  dev.off()

  # Computing the mean Fst for each pop and SDs
  npop <- length(unique(within$FID))
  pop <- unique(within$FID)
  mean_pop <- array(NA, length(pop))
  sd_pop <- array(NA, length(pop))
  for (i in 1:length(pop)) {
    a <- subset(pwFst, pwFst$V1 == pop[i])
    b <- subset(pwFst, pwFst$V2 == pop[i])
    c <- rbind.data.frame(a, b)
    mean_pop[i] <- mean(c$V3)
    sd_pop[i] <- sd(c$V3)
    }
  mean_pop_df <- data.frame(Pop = pop, Mean_Fst = mean_pop, SD = sd_pop, X = NA, Y = NA)
  
  # Adding the coordinates
  coordInd <- as.data.frame(data.table::fread(coordInd, header = TRUE))
  colnames(coordInd) <- c("id", "long", "lat")
  X <- by(data = coordInd$long, INDICES = coordInd$id, FUN = mean)
  pop_coord <- names(X)
  X <- as.vector(X)
  Y <- as.vector(by(data = coordInd$lat, INDICES = coordInd$id, FUN = mean))
  pop_coord_df <- data.frame(Pop = pop_coord, X = X, Y = Y)
  
  for (i in 1:nrow(mean_pop_df)) {
    mean_pop_df$X[i] <- pop_coord_df$X[which(mean_pop_df$Pop[i] == pop_coord_df$Pop)]
    mean_pop_df$Y[i] <- pop_coord_df$Y[which(mean_pop_df$Pop[i] == pop_coord_df$Pop)]
  }
  
  return(mean_pop_df)

}