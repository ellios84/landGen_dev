geneticDiversityMap <- function(diversity_df = NULL, 
                                diversity_stat = NULL, 
                                coordInd = NULL, 
                                elev = "wc2.1_30s_elev.tif", 
                                extent_factor = NULL, 
                                scalebar_pos = NULL) {
  
  # Loading individual coordinates
  coordInd <- as.data.frame(data.table::fread(coordInd, header = TRUE))
  colnames(coordInd) <- c("id", "long", "lat")
  
  X <- by(data = coordInd$long, INDICES = coordInd$id, FUN = mean)
  pop_coord <- names(X)
  X <- as.vector(X)
  Y <- as.vector(by(data = coordInd$lat, INDICES = coordInd$id, FUN = mean))
  pop_coord_df <- data.frame(Pop = pop_coord, X = X, Y = Y)
  
  # Loading elevation and deriving the map background
  elev <- raster::raster(elev)
  ext <- raster::extent(range(pop_coord_df$X)[1] - extent_factor, range(pop_coord_df$X)[2] + extent_factor, range(pop_coord_df$Y)[1] - extent_factor, range(pop_coord_df$Y)[2] + extent_factor)
  elev <- raster::crop(x = elev, y = ext)
  slope <- raster::terrain(x = elev, opt = "slope")
  aspect <- raster::terrain(x = elev, opt = "aspect")
  hill <- raster::hillShade(slope = slope, aspect = aspect)
  col_elev <- gray.colors(n = 256, start = 0.2, end = .99, alpha = 1, rev = TRUE, gamma = 2.5)
  col_hill <- gray.colors(n = 256, start = 0.2, end = .99, alpha = .2, rev = FALSE, gamma = 2.5)
  colPalAlt <- gray.colors(n = 256, start = 0.2, end = .99, alpha = 1, rev = TRUE, gamma = 2.5)
  colPalHill <- gray.colors(n = 256, start = 0.2, end = .99, alpha = .2, rev = FALSE, gamma = 2.5)
  
  # Adding the population coordinates
  diversity_df$X <- NA
  diversity_df$Y <- NA
  for(i in 1:nrow(pop_coord_df)) {
    diversity_df$X[which(pop_coord_df$Pop[i] == diversity_df$pop)] <- pop_coord_df$X[i]
    diversity_df$Y[which(pop_coord_df$Pop[i] == diversity_df$pop)] <- pop_coord_df$Y[i]
  }
  
  diversityQuantile <- quantile(diversity_df[, diversity_stat])
  diversity_df$Col <- NA
  diversity_df$Col[which(diversity_df[, diversity_stat] >= diversityQuantile[1] & diversity_df[, diversity_stat] < diversityQuantile[2])] <- "#007FFF"
  diversity_df$Col[which(diversity_df[, diversity_stat] >= diversityQuantile[2] & diversity_df[, diversity_stat] < diversityQuantile[3])] <- "#99EDFF"
  diversity_df$Col[which(diversity_df[, diversity_stat] >= diversityQuantile[3] & diversity_df[, diversity_stat] < diversityQuantile[4])] <- "#FFC34C"
  diversity_df$Col[which(diversity_df[, diversity_stat] >= diversityQuantile[4])] <- "#FF7F00"
  
  SdQuantile <- quantile(diversity_df[, paste0(diversity_stat, "SD")])
  diversity_df$SD_cex <- NA
  diversity_df$SD_cex[which(diversity_df[, paste0(diversity_stat, "SD")] >= SdQuantile[1] & diversity_df[, paste0(diversity_stat, "SD")] < SdQuantile[2])] <- 1
  diversity_df$SD_cex[which(diversity_df[, paste0(diversity_stat, "SD")] >= SdQuantile[2] & diversity_df[, paste0(diversity_stat, "SD")] < SdQuantile[3])] <- 1.4
  diversity_df$SD_cex[which(diversity_df[, paste0(diversity_stat, "SD")] >= SdQuantile[3] & diversity_df[, paste0(diversity_stat, "SD")] < SdQuantile[4])] <- 1.8
  diversity_df$SD_cex[which(diversity_df[, paste0(diversity_stat, "SD")] >= SdQuantile[4])] <- 2.2
  
  jpeg(paste0(diversity_stat, "-map.jpeg"), width = 10, height = 8, units = 'in', res = 800)
  layout(mat = matrix(c(1, 2), 1), widths = c(1, 0.2), heights = c(1, 1))
  par(mar = c(0, 0, 0, 0), oma = c(4.5, 4.5, 4.5, 4.5))
  
  raster::image(elev, col = colPalAlt, xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
  raster::image(hill, col = colPalHill, xaxt = "n", yaxt = "n", bty="n", xlab = "", ylab = "", add = TRUE)
  points(diversity_df$X, diversity_df$Y, pch = 21, bg = diversity_df$Col, cex = diversity_df$SD_cex)
  prettymapr::addscalebar(style = "bar", plotepsg ="4326", plotunit = "Km", pos = "bottomright", padin = c(0.1, 0.1), label.cex = 1, widthhint = 0.2, label.col = "black", bar.cols = c("black", "black"), htin = 0.02, linecol = "black")
  box()
  
  plot(1:10, 1:10, axes = FALSE, t = "n")
  if (diversity_stat == "FIS") {
    text(5.9, 8.6, expression(italic(F)[IS]), cex = 1.3)
    }
  else {
    text(5.9, 8.6, diversity_stat, cex = 1.3)
  }
  
  points(3, 8, pch = 22, cex = 4, bg = "#FF7F00")
  points(3, 7.6, pch = 22, cex = 4, bg = "#FFC34C")
  points(3, 7.2, pch = 22, cex = 4, bg = "#99EDFF")
  points(3, 6.8, pch = 22, cex = 4, bg = "#007FFF")
  text(6, 6.6, sprintf("%5.3f", round(diversityQuantile[1], digits = 3)))
  text(6, 7, sprintf("%5.3f", round(diversityQuantile[2], digits = 3)))
  text(6, 7.4, sprintf("%5.3f", round(diversityQuantile[3], digits = 3)))
  text(6, 7.8, sprintf("%5.3f", round(diversityQuantile[4], digits = 3)))
  text(6, 8.2, sprintf("%5.3f", round(diversityQuantile[5], digits = 3)))
  
  text(4, 5.8, "Std dev", cex = 1.3)
  points(3, 5.2, pch = 22, cex = 4)
  points(3, 5.2, cex = 2.2)
  points(3, 4.78, pch = 22, cex = 4)
  points(3, 4.78, cex = 1.8)
  points(3, 4.36, pch = 22, cex = 4)
  points(3, 4.36, cex = 1.4)
  points(3, 3.94, pch = 22, cex = 4)
  points(3, 3.94, cex = 1)
  text(6, 3.74, sprintf("%5.3f", round(SdQuantile[1], digits = 3)))
  text(6, 4.16, sprintf("%5.3f", round(SdQuantile[2], digits = 3)))
  text(6, 4.58, sprintf("%5.3f", round(SdQuantile[3], digits = 3)))
  text(6, 5, sprintf("%5.3f", round(SdQuantile[4], digits = 3)))
  text(6, 5.4, sprintf("%5.3f", round(SdQuantile[5], digits = 3)))
  
  dev.off()

}