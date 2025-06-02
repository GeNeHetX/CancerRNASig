#' Title applySignature
#' 
#' @param newexp gene expression matrix/dataframe, sample in columns, gene in rows
#' @param geneSymbols gene symbols, a vector of same length as the number of rows in newex
#' @param signature signature to apply, can be Gempred, Puleo, Purist or McpCount
#' @param norm normlisation to apply if needed, can be none(by default), uq or vst 
#' @param scale scale, can be raw, sc=sample center, gc=gene center, gsc=gene scale center, ssc=sample scale center (none by default)
#'
#' @return a data frame with row names as colnames of newexp, columns are different according the signatures, can be signature score or/and signature conclusion
#'
#' @examples 
#' applySignature(raw_counts,geneannot,signature="Gempred",norm="uq",scale="ssc")
#' @export


applySignature = function (matrix, geneSymbols, signature = c("Gempred", "Puleo", "Purist", "McpCount"), norm="none", scale="none"){
  
  # Input validation
  if (!is.matrix(matrix)){
    matrix = as.matrix(matrix)
  }
  if (!is.character(geneSymbols)) stop("geneSymbols must be a character vector.")
  sig <- match.arg(signature)  # Verify that the signature argument is valid
  
  # Apply normalisation
  if(norm=="uq"){
    print("Launching UQ normalisation")
    data = qutils::UQnorm(matrix)
  }
  else if(norm=="vst"){
    print("Launching VST normalisation")
    data = .VSTnorm(matrix)
  }
  else{
    print("No normalisation applied")
    data = matrix
  }
  
  # Apply scale
  data = CancerRNASig:::.qNormalize(data, scale)
  
  # Apply signature
  if (sig=="McpCount"){
    print("Launching McpCount prediction")
    data <- getUniqueGeneMat(data, geneSymbols, rowMeans(data))
    res = CancerRNASig::MCPcounter(data)
  }
  else if (sig=="Puleo"){
    print("Launching Puleo prediction")
    res = CancerRNASig::qProjICA(getUniqueGeneMat(data, geneSymbols, rowMeans(data)))
  }
  else if (sig=="Purist"){
    print("Launching Purist prediction")
    res = CancerRNASig::purist(data,geneSymbols)
    #names(res) <- c("PuristScore","PuristPrediction")
  }
  else if (sig=="Gempred"){
    print("Launching GemPred 2 prediction")
    res = CancerRNASig::GemPred(data, geneSymbols)
  }
  
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