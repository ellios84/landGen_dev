pairwiseFstMap <- function(MeanPairwiseFst = NULL, elev = "wc2.1_30s_elev.tif", extent_factor = NULL, scalebar_pos = NULL) {
  
  elev <- raster::raster(elev)
  ext <- raster::extent(range(MeanPairwiseFst$X)[1] - extent_factor, range(MeanPairwiseFst$X)[2] + extent_factor, 
                        range(MeanPairwiseFst$Y)[1]- extent_factor, range(MeanPairwiseFst$Y)[2] + extent_factor)
  elev <- raster::crop(x = elev, y = ext)
  slope <- raster::terrain(x = elev, opt = "slope")
  aspect <- raster::terrain(x = elev, opt = "aspect")
  hill <- raster::hillShade(slope = slope, aspect = aspect)
  col_elev <- gray.colors(n = 256, start = 0.2, end = .99, alpha = 1, rev = TRUE, gamma = 2.5)
  col_hill <- gray.colors(n = 256, start = 0.2, end = .99, alpha = .2, rev = FALSE, gamma = 2.5)
  colPalAlt <- gray.colors(n = 256, start = 0.2, end = .99, alpha = 1, rev = TRUE, gamma = 2.5)
  colPalHill <- gray.colors(n = 256, start = 0.2, end = .99, alpha = .2, rev = FALSE, gamma = 2.5)
  
  FstQuantile <- quantile(MeanPairwiseFst$Mean_Fst)
  MeanPairwiseFst$Col <- NA
  MeanPairwiseFst$Col[which(MeanPairwiseFst$Mean_Fst >= FstQuantile[1] & MeanPairwiseFst$Mean_Fst < FstQuantile[2])] <- "#007FFF"
  MeanPairwiseFst$Col[which(MeanPairwiseFst$Mean_Fst >= FstQuantile[2] & MeanPairwiseFst$Mean_Fst < FstQuantile[3])] <- "#99EDFF"
  MeanPairwiseFst$Col[which(MeanPairwiseFst$Mean_Fst >= FstQuantile[3] & MeanPairwiseFst$Mean_Fst < FstQuantile[4])] <- "#FFC34C"
  MeanPairwiseFst$Col[which(MeanPairwiseFst$Mean_Fst >= FstQuantile[4])] <- "#FF7F00"
  
  SdQuantile <- quantile(MeanPairwiseFst$SD)
  MeanPairwiseFst$SD_cex <- NA
  MeanPairwiseFst$SD_cex[which(MeanPairwiseFst$SD >= SdQuantile[1] & MeanPairwiseFst$SD < SdQuantile[2])] <- 1
  MeanPairwiseFst$SD_cex[which(MeanPairwiseFst$SD >= SdQuantile[2] & MeanPairwiseFst$SD < SdQuantile[3])] <- 1.4
  MeanPairwiseFst$SD_cex[which(MeanPairwiseFst$SD >= SdQuantile[3] & MeanPairwiseFst$SD < SdQuantile[4])] <- 1.8
  MeanPairwiseFst$SD_cex[which(MeanPairwiseFst$SD >= SdQuantile[4])] <- 2.2
  
  jpeg("pairwiseFst-map.jpeg", width = 10, height = 8, units = 'in', res = 800)
  layout(mat = matrix(c(1, 2), 1), widths = c(1, 0.2), heights = c(1, 1))
  par(mar = c(0, 0, 0, 0), oma = c(4.5, 4.5, 4.5, 4.5))
  
  raster::image(elev, col = colPalAlt, xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
  raster::image(hill, col = colPalHill, xaxt = "n", yaxt = "n", bty="n", xlab = "", ylab = "", add = TRUE)
  points(MeanPairwiseFst$X, MeanPairwiseFst$Y, pch = 21, bg = MeanPairwiseFst$Col, cex = MeanPairwiseFst$SD_cex)
  prettymapr::addscalebar(style = "bar", plotepsg ="4326", plotunit = "Km", pos = "bottomright", padin = c(0.1, 0.1), label.cex = 1, widthhint = 0.2, label.col = "black", bar.cols = c("black", "black"), htin = 0.02, linecol = "black")
  box()
  
  plot(1:10, 1:10, axes = FALSE, t = "n")
  text(5.9, 8.6, expression(paste("Pairwise ", italic(F)[ST])), cex = 1.3)
  points(3, 8, pch = 22, cex = 4, bg = "#FF7F00")
  points(3, 7.6, pch = 22, cex = 4, bg = "#FFC34C")
  points(3, 7.2, pch = 22, cex = 4, bg = "#99EDFF")
  points(3, 6.8, pch = 22, cex = 4, bg = "#007FFF")
  text(6, 6.6, sprintf("%5.3f", round(FstQuantile[1], digits = 3)))
  text(6, 7, sprintf("%5.3f", round(FstQuantile[2], digits = 3)))
  text(6, 7.4, sprintf("%5.3f", round(FstQuantile[3], digits = 3)))
  text(6, 7.8, sprintf("%5.3f", round(FstQuantile[4], digits = 3)))
  text(6, 8.2, sprintf("%5.3f", round(FstQuantile[5], digits = 3)))
  
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