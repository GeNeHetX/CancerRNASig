#' Title callSignature
#' 
#' @param newexp gene expression matrix/dataframe, sample in columns, gene in rows
#' @param geneSymbols gene symbols, a vector of same length as the number of rows in newex
#' @param signature signature to apply, can be Gempred, GempredBiopsy, Puleo, McpCount or Purist
#' @param toNorm normlisation to apply if needed, can be none (by default), uq or vst 
#' @param toScale scale to apply if needed, can be none (by default), sc=sample center, gc=gene center, gsc=gene scale center, ssc=sample scale center
#'
#' @return a data frame with row names as colnames of newexp, columns are different according the signatures, can be signature score or/and signature conclusion
#'
#' @examples 
#' callSignature(raw_counts,geneannot,signature="Gempred",toNorm="uq",scale="ssc")
#' @export


callSignature = function (matrix, geneSymbols, signature = c("Gempred", "GempredBiopsy", "Puleo", "McpCount", "Purist"), toNorm="none", toScale="none"){
  
  # -- Input validation & coercions -------------------------------------------------
  if (!is.matrix(matrix))
    matrix <- as.matrix(matrix)
  if (all(grepl("^ENSG", rownames(data))))
    stop("Row names of input matrix should be ENSG ids.")
  if (!is.character(geneSymbols))
    stop("geneSymbols must be a character vector.")

  sig     <- match.arg(signature)
  toNorm  <- match.arg(toNorm)
  toScale <- match.arg(toScale)
  
  # -- Normalisation ----------------------------------------------------------
  data <- switch(toNorm,
    uq   = { message("Launching UQ normalisation") ; qutils::UQnorm(matrix) },
    vst  = { message("Launching VST normalisation"); vst(matrix)       },
    none = { message("No normalisation applied")   ; matrix                 }
  )

  # -- Scaling -----------------------------------------------------------------
  data <- .qNormalize(data, toScale)

  # -- GeneSymbol -----------------------------------------------------------------

  # -- Signatures ---------------------------------------------------
  res <- switch(sig,
    McpCount = {
      message("Launching McpCount prediction")
      #data_u <- getUniqueGeneMat(data, geneSymbols, rowMeans(data))
      .mcpcount(data, geneSymbols)
    },
    Puleo = {
      message("Launching Puleo prediction")
      #raw data, no normalisation
      .qProjICA(getUniqueGeneMat(data, geneSymbols, rowMeans(data)))
    },
    Purist = {
      message("Launching Purist prediction")
      .purist(data, geneSymbols)
    },
    ## add estimate signature 

    #need ENSG or Gene Symbols
    Gempred = {
      message("Launching GemPred prediction")
      .gemPred(data, geneSymbols)
    },
    ## add GempredBiopsy - need ENSG
    GempredBiopsy = {
      message("Launching GemPred Biopsy prediction")
      .GemPred_biopsy(data, gnt="raw",pnt="raw")
    },
    
    stop("Unreachable state – unknown signature. Please choose from 'Gempred', 'GempredBiopsy', 'Puleo', 'McpCount' or 'Purist'.")  # sécurité
  )
  
  res
}