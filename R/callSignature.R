#' Title callSignature
#' 
#' @param newexp gene expression matrix/dataframe, sample in columns, gene in rows
#' @param geneSymbols gene symbols, a vector of same length as the number of rows in newex
#' @param signature signature to apply, can be Gempred, Puleo, McpCount or Purist
#' @param toNorm normlisation to apply if needed, can be none(by default), uq or vst 
#' @param toScale scale, can be raw, sc=sample center, gc=gene center, gsc=gene scale center, ssc=sample scale center (none by default)
#'
#' @return a data frame with row names as colnames of newexp, columns are different according the signatures, can be signature score or/and signature conclusion
#'
#' @examples 
#' callSignature(raw_counts,geneannot,signature="Gempred",toNorm="uq",scale="ssc")
#' @export


callSignature = function (matrix, geneSymbols, signature = c("Gempred", "Puleo", "McpCount", "Purist"), toNorm="none", toScale="none"){
  
  # -- Input validation & coercions -------------------------------------------------
  if (!is.matrix(matrix))
    matrix <- as.matrix(matrix)
  if (!is.character(geneSymbols))
    stop("geneSymbols must be a character vector.")

  sig     <- match.arg(signature)
  toNorm  <- match.arg(toNorm)
  toScale <- match.arg(toScale)
  
  # -- Normalisation ----------------------------------------------------------
  data <- switch(toNorm,
    uq   = { message("Launching UQ normalisation") ; qutils::UQnorm(matrix) },
    vst  = { message("Launching VST normalisation"); .VSTnorm(matrix)       },
    none = { message("No normalisation applied")   ; matrix                 }
  )

  # -- Scaling -----------------------------------------------------------------
  data <- .qNormalize(data, toScale)

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
    Gempred = {
      message("Launching GemPred prediction")
      .gemPred(data, geneSymbols)
    },
    ## add gempred_biopsi
    
    stop("Unreachable state – unknown signature. Please choose from 'Gempred', 'Puleo', 'McpCount' or 'Purist'.")  # sécurité
  )
  
  res
}

## ADD to qutils package like UQnorm function ?
.VSTnorm = function(counts){
    conditions <- character(0)
    mid = as.integer(dim(counts)[2]/2)
    conditions <- c(rep("A",mid),rep("B",dim(counts)[2]-mid))
    metadata <- data.frame(condition = factor(conditions))
    rownames(metadata) <- colnames(counts)
    dds <- DESeqDataSetFromMatrix(
      countData = counts,
      colData = metadata,
      design = ~ condition
    )
    dds <- DESeq(dds)
    dds =  vst(dds)
    vst = assay(dds)
    return(vst)
  }