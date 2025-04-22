library(XML)
library(curl)
library(dplyr)
library(biomaRt)
library(xml2)

# Parameters for gene annotation
ensembl = useEnsembl(biomart = "ensembl")
ensembl = useMart("ensembl", dataset="oaries_gene_ensembl",
                  ensemblRedirect = FALSE)


sign <- as.data.frame(as.character(unique(outliers$SNP)))
colnames(sign) <- "SNPid"

map <- read.table("new_IMAGE_sheepHD_chrok_ok.map", stringsAsFactors = F)

match_map <- function(x) {
  t <- which(as.character(sign$SNPid[x]) == map$V2)
  t <- map[t, c(1, 4)]
  return(t)
}

res.ann <- plyr::adply(1:nrow(sign), 1, match_map)
res.ann <- res.ann[, 2:3]
res.ann <- cbind.data.frame(sign, res.ann)
colnames(res.ann) <- c("Marker", "Chr", "Pos")

gw <- 20000
res.ann$Start <- res.ann$Pos - gw
res.ann$End <- res.ann$Pos + gw

res.ann <- res.ann[order(res.ann$Chr), ]

anf <- function(x) {
  anfx <- biomaRt::getBM(attributes = c("ensembl_gene_id", 
                                        "wikigene_name",
                                        "gene_biotype",
                                        "external_gene_name", 
                                        "external_gene_source",
                                        "start_position", "end_position", 
                                        # "minor_allele",
                                        "description", 
                                        "go_id", 
                                        "goslim_goa_accession",
                                        "goslim_goa_description"),
                         filters = c("chromosome_name","start","end"),
                         values = list(res.ann$Chr[x], res.ann$Start[x], res.ann$End[x]),
                         mart = ensembl)
  if(dim(anfx)[1] != 0){
    anfx <- cbind.data.frame(res.ann[x, ], anfx, row.names = NULL)
    return(anfx)
  }
}

res.ann <- plyr::adply(c(1:nrow(res.ann)), 1, anf)
