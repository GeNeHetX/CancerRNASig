#' Title GemPred_Biopsy
#'
#' @description WARNING: This tool is intended for research purposes only. It is not intended for diagnostic use or clinical decision-making.
#' 
#' @param newexp gene expression matrix/dataframe, sample in columns, gene in rows
#'
#' @return a data frame with row names as colnames of newexp, first column is gempred score
#' @keywords internal
#' 
#' @examples
#' # Example of how to call GemPred_Biopsy
#' # result <- CancerRNASig:::.GemPred_Biopsy(newexp = gene_expression_matrix)
#'

.GemPred_biopsy = function(newexp,  gnt = "sc", pnt = "raw") {
  
  data(t_GemPred_Biopsy)
  sig <- as.vector(t_GemPred_Biopsy$weight) ; names(sig) <- t_GemPred_Biopsy$ENG_ID
  proj <- as.data.frame(CancerRNASig:::.qProjICA(newexp, ICAgw = as.matrix(sig), geneNormType = gnt, projNormType = pnt, ming = 1))
  colnames(proj) <- "newsig"
  
  return(proj)
  
}