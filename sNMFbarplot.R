# Important note: the 'sNMFbarplot' function requires as input the object returned by plotCedecay()
# snmf_i = sNMF results
# plotK = c(..., ...) --> the numebr of Ks you want to plot
# coordInd = as for pairwiseFst()

sNMFbarplot <- function(snmf_res = NULL, CEdecay = NULL, plotK = NULL, 
                        grplab = NULL, sortind = NA,  
                        grplabsize = NULL, grplabangle = NULL, 
                        grplabspacer = NULL, grplabheight = NULL, grplabpos = NULL,   
                        plot_height = 25, plot_width = 30, 
                        grpmean = FALSE, barbordersize = NULL, linesize = NULL, 
                        legendkeysize = NULL, legendtextsize = NULL, legendmargin = c(0.5, 0.5, 0.5, 0),
                        ticksize = NULL, splabsize = NULL, pointsize = NULL, divsize = NULL, panelspacer = NULL) {
  
  CEdecay <- CEdecay[[1]]
  bestruns <- as.vector(by(data = CEdecay$CE, INDICES = CEdecay$K, FUN = which.min))
  names(bestruns) <- unique(CEdecay$K)
  bestruns <- bestruns[plotK]
  
  q_ls <- list()
  for (i in 1:length(bestruns)) {
    q_i <- LEA::Q(snmf_res, K = as.numeric(names(bestruns)[i]), run = bestruns[i]) 
    colnames(q_i) <- paste0("Cluster", 1:(as.numeric(names(bestruns)[i])))
    q_ls[[i]] <- as.data.frame(q_i)
  }

  # Checking the right format
  q_ls <- pophelper::as.qlist(q_ls)
  print(pophelper::summariseQ(pophelper::tabulateQ(q_ls)))
  
  # Aligning Ks (likewise Clump)
  q_ls <- pophelper::alignK(qlist = q_ls, type = "across")

  # Labelling
  grplab <- read.table(file = grplab, header = TRUE)

  # Cluster colors
  set.seed(123)
  mycol <- sample(x = colorBlindness::PairedColor12Steps, size = max(plotK), replace = FALSE)
  
  p1 <- pophelper::plotQ(q_ls, imgoutput = "join", returnplot = TRUE, exportplot = FALSE, basesize = 11,
                         grplab = grplab, grplabangle = grplabangle, grplabsize = grplabsize, 
                         grplabspacer = grplabspacer, grplabheight = grplabheight, grplabpos = grplabpos,
                         barbordersize = barbordersize, barbordercolour = "white",
                         linecol = "black", linealpha = 1, 
                         grpmean = grpmean, pointsize = pointsize, 
                         linesize = linesize, panelspacer = panelspacer,
                         divsize = divsize, divtype = 1, splab = paste0("K = ", plotK),  
                         sharedindlab = FALSE, clustercol = mycol, sortind = sortind, 
                         splabsize = splabsize, showyaxis = FALSE, showticks = TRUE, ticksize = ticksize, 
                         showindlab = FALSE,
                         showlegend = TRUE, legendpos = "left", legendkeysize = legendkeysize, 
                         legendtextsize = legendtextsize, legendrow = 1, legendspacing = 30, 
                         legendmargin = legendmargin)
  
  ggplot2::ggsave(filename = "sNMFbarplot.jpeg", limitsize = FALSE, 
                  plot = gridExtra::grid.arrange(p1$plot[[1]]),
                  dpi = 320, height = plot_height, width = plot_width, device = "jpeg")
}
