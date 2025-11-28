#' Title MCPcounter
#'
#' @param newexp gene expression matrix/dataframe, sample in columns, gene in rows
#' @param geneSymbols gene symbols, a vector of same length as the number of rows in newex
#'
#' @return a data frame with samples in row and MCPcounter population quantification in columns
#' @keywords internal
#'
#' @examples
.mcpcount <- function(newexp, geneSymbols) {
  if (nrow(newexp) != length(geneSymbols)) {
    stop("geneSymbols should be a vector of gene symbols exactly corresponding to each row of the newexp dataset")
  }
  # genemat=CancerRNASig:::.getugm(newexp,geneSymbols,matrixStats::rowSds(as.matrix(newexp),useNames=T))
  genemat <- newexp

  markers.names <- c(
    "Tcells", "CD8Tcells", "Cytotox.lymph",
    "NK", "B.lineage", "Mono.lineage", "Myeloid.dendritic",
    "Neutrophils", "Endothelial", "Fibroblasts"
  )

  markers <- CancerRNASig:::mcpgenes
  features <- subset(markers, markers$HUGO.symbols %in% rownames(genemat))
  features <- split(features[, "HUGO.symbols"], features[, "Cell.population"])
  missing.populations <- setdiff(markers.names, names(features))
  features <- features[intersect(markers.names, names(features))]
  if (length(missing.populations) > 0) {
    warning(paste(
      "Found no markers for population(s):",
      paste(missing.populations, collapse = ", ")
    ))
  }

  res <- as.data.frame(do.call(cbind, lapply(features, function(x) {
    apply(genemat[intersect(row.names(genemat), x), , drop = F], 2,
      mean,
      na.rm = T
    )
  })))
  res
}
