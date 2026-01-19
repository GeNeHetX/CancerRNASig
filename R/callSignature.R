#' Title callSignature
#'
#' @param matrix gene expression matrix/dataframe, sample in columns, gene in rows
#' @param geneSymbols gene symbols, a vector of same length as the number of rows in matrix
#' @param signature signature to apply, can be Gempred, tGempred, Puleo, pdacMolGrad, Mcpcount or Purist
#' @param normType normlisation to apply if needed, can be raw (by default), uq or vst
#' @param scaleType scale to apply if needed, can be raw (by default), sc=sample center, gc=gene center, gsc=gene scale center, ssc=sample scale center
#' @param ... Additional arguments.
#'
#' @details
#' The \code{...} argument allows passing additional parameters to specific
#' signature implementations.
#'
#' Currently supported extra parameters include:
#' \itemize{
#'   \item \strong{mingenes}: Minimum number of genes in common required to compute
#'   the \code{"Gempred"} signature (default: 500).
#' }
#'
#' Additional parameters may be ignored depending on the selected signature.
#'
#' @return a data frame with row names as colnames of newexp, columns are different according the signatures, can be signature score or/and signature conclusion
#'
#' @examples
#' callSignature(matrix, geneannot, signature = "Gempred", normType = "uq", scale = "ssc", ...)
#' @export


callSignature <- function(matrix,
                          geneSymbols = NULL,
                          signature = NULL,
                          normType = "raw",
                          scaleType = "raw",
                          ...) {
  # -- Input validation & coercions -------------------------------------------------
  if (!is.matrix(matrix)) {
    matrix <- as.matrix(matrix)
  }

  if (!is.character(geneSymbols) && !is.null(geneSymbols)) {
    stop("geneSymbols must be a character vector.")
  }

  if (is.null(geneSymbols)) {
    if (all(grepl("^ENSG", rownames(matrix))) && !(signature %in% c("Gempred", "tGempred"))) {
      stop("Matrix rownames detected as Ensembl Gene IDs. Whitout geneSymbols provided, only gempred and tGempred signatures can be applied.")
    } else {
      geneSymbols <- rownames(matrix)
    }
  }

  if (is.null(signature)) stop("You must specify a signature: 'Gempred', 'tGempred', 'Puleo', 'pdacMolGrad', 'Mcpcount', or 'Purist'.")
  sig <- match.arg(signature, choices = c("Gempred", "tGempred", "Puleo", "Mcpcount", "Purist", "pdacMolGrad"))

  normType <- match.arg(normType, choices = c("raw", "uq", "vst"))
  scaleType <- match.arg(scaleType, choices = c("raw", "gc", "sc", "gsc", "ssc"))

  # -- Normalisation ----------------------------------------------------------
  data <- switch(normType,
    uq = {
      message("Launching UQ normalisation")
      CancerRNASig:::.UQnorm(matrix)
    },
    vst = {
      message("Launching VST normalisation")
      DESeq2::vst(matrix)
    },
    raw = {
      message("No normalisation applied")
      matrix
    }
  )

  # -- Scaling -----------------------------------------------------------------
  data <- CancerRNASig:::.qNormalize(data, scaleType)

  # -- Unique Gene Matrix -----------------------------------------------------
  # data_unique <- qutils::getUniqueGeneMat(data, geneSymbols, rowMeans(data))
  if (!is.null(geneSymbols)) {
    data_unique <- CancerRNASig:::genesymUniqExp(data, geneSymbols, rowMeans)
  } else {
    data_unique <- data
  }


  # -- Signatures ---------------------------------------------------
  res <- switch(sig,
    Mcpcount = {
      message("Launching Mcpcount prediction")
      CancerRNASig::MCPcounter(data_unique, rownames(data_unique))
    },
    Puleo = {
      message("Launching Puleo prediction")
      # raw data, no normalisation
      CancerRNASig::qProjICA(data_unique, ...)
    },
    Purist = {
      message("Launching Purist prediction")
      CancerRNASig::purist(data_unique, rownames(data_unique))
    },
    Gempred = {
      message("Launching Gemred prediction")
      CancerRNASig::gemPred(data_unique, geneSymbols, ...)
    },
    tGempred = {
      message("Launching tGempred prediction")
      CancerRNASig::GemPred_biopsy(data_unique, gnt = "raw", pnt = "raw")
    },
    pdacMolGrad = {
      message("Launching pdacMolGrad")
      CancerRNASig::projectMolGrad(data_unique, row.names(data_unique), "raw")
    },

    ## add estimate signature
    stop("Unreachable state – unknown signature. Please choose from 'Gempred', 'tGempred', 'Puleo', 'pdacMolGrad','Mcpcount' or 'Purist'.") # sécurité
  )

  res
}
