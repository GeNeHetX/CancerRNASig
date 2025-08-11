#' Title purist
#'
#' @param newexp gene expression matrix/dataframe, sample in columns, gene in rows
#' @param geneSymbols gene symbols, a vector of same length as the number of rows in newex
#'
#' @return a data frame with row names as colnames of newexp, first column is purist subtype, second is purist score
#' @keywords internal
#'
#' @examples
#' 

.purist=function(newexp,geneSymbols){

  if(nrow(newexp)!= length(geneSymbols)){
    stop("geneSymbols should be a vector of gene symbols exactly corresponding to each row of the newexp dataset")
  }
  expg=qutils::getUniqueGeneMat(newexp,geneSymbols,matrixStats::rowSds(as.matrix(newexp),useNames = TRUE))

  inter=-6.815


  puristcoef=data.frame(
    a=c("GPR87","KRT6A","BCAR3","PTGES","ITGA3","C16orf74","S100A2","KRT5"),
    b=c("REG4","ANXA10","GATA6","CLDN18","LGALS4","DDC","SLC40A1","CLRN3"),
    w=c(1.994,2.031,1.618,0.922,1.059,0.929,2.505,0.485),stringsAsFactors = F)
  puristg=c(puristcoef$a,puristcoef$b)


  if(!all(puristg %in% rownames(expg))){
    stop(paste("Missing some genes:",paste(puristg[which(!puristg %in% rownames(expg))],collapse=", ")))
  }



  prednum=setNames(rowSums(apply(puristcoef,1,function(p){(expg[p[1],]>expg[p[2],])*as.numeric(p[3])}))+inter,colnames(expg))
  preds=setNames(factor(c("classic","basal")[ (prednum>0)+1]),colnames(expg))


  data.frame(NumPurist=prednum,
             Purist=preds,
             row.names=colnames(expg),stringsAsFactors = F)
}
