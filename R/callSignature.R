#' Title callSignature
#' 
#' @param matrix gene expression matrix/dataframe, sample in columns, gene in rows
#' @param geneSymbols gene symbols, a vector of same length as the number of rows in matrix
#' @param signature signature to apply, can be Gempred, tGempred, Puleo, Mcpcount or Purist
#' @param normType normlisation to apply if needed, can be raw (by default), uq or vst 
#' @param scaleType scale to apply if needed, can be raw (by default), sc=sample center, gc=gene center, gsc=gene scale center, ssc=sample scale center
#'
#' @return a data frame with row names as colnames of newexp, columns are different according the signatures, can be signature score or/and signature conclusion
#'
#' @examples 
#' callSignature(matrix,geneannot,signature="Gempred",normType="uq",scale="ssc")
#' @export


callSignature = function (matrix, geneSymbols, signature = NULL, normType="raw", scaleType="raw"){
  
  # -- Input validation & coercions -------------------------------------------------
  if (!is.matrix(matrix))
    matrix <- as.matrix(matrix)
  if (!all(grepl("^ENSG", rownames(matrix))))
    stop("Row names of input matrix should be ENSG ids.")
  if (!is.character(geneSymbols))
    stop("geneSymbols must be a character vector.")
  
  if (is.null(signature)) stop("You must specify a signature: 'Gempred', 'tGempred', 'Puleo', 'Mcpcount', or 'Purist'.")
  sig <- match.arg(signature, choices = c("Gempred", "tGempred", "Puleo", "Mcpcount", "Purist"))
  
  normType  <- match.arg(normType, choices = c("raw", "uq", "vst"))
  scaleType <- match.arg(scaleType, choices = c("raw", "gc", "sc", "gsc", "ssc"))
  
  # -- Normalisation ----------------------------------------------------------
  data <- switch(normType,
    uq   = { message("Launching UQ normalisation") ; qutils::UQnorm(matrix) },
    vst  = { message("Launching VST normalisation"); DESeq2::vst(matrix)       },
    raw = { message("No normalisation applied")   ; matrix                 }
  )

  # -- Scaling -----------------------------------------------------------------
  data <- qutils::qNormalize(data, scaleType)

  # -- Signatures ---------------------------------------------------
  res <- switch(sig,
    Mcpcount = {
      message("Launching Mcpcount prediction")
      #data_u <- getUniqueGeneMat(data, geneSymbols, rowMeans(data))
      CancerRNASig:::.mcpcount(data, geneSymbols)
    },
    Puleo = {
      message("Launching Puleo prediction")
      #raw data, no normalisation
      .qProjICA(qutils::getUniqueGeneMat(data, geneSymbols, rowMeans(data)))
    },
    Purist = {
      message("Launching Purist prediction")
      CancerRNASig:::.purist(data, geneSymbols)
    },
    ## add estimate signature 

    #need ENSG or Gene Symbols
    Gempred = {
      message("Launching Gemred prediction")
      CancerRNASig:::.gemPred(data, geneSymbols)
    },
    ## add GempredBiopsy - need ENSG
    tGempred = {
      message("Launching tGempred prediction")
      CancerRNASig:::.GemPred_biopsy(data, gnt="raw",pnt="raw")
    },
    
    stop("Unreachable state – unknown signature. Please choose from 'Gempred', 'tGempred', 'Puleo', 'Mcpcount' or 'Purist'.")  # sécurité
  )
  
  res
}