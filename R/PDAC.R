
#' Title Purist
#'
#' @param newexp gene expression matrix/dataframe, sample in columns, gene in rows.
#'  Gene symbols are required
#'
#' @return a data frame with row names as colnames of newexp, first column is purist subtype, second is purist score
#' @export
#'
#' @examples
# PDAC=function(newexp,system=c("puleo","pamg","purist")){

#   .checkGeneExp(newexp)


# }


# #' Title Purist
# #'
# #' @param newexp gene expression matrix/dataframe, sample in columns, gene in rows.
# #'  Gene symbols are required
# #'
# #' @return a data frame with row names as colnames of newexp, first column is purist subtype, second is purist score
# #' @export
# #'
# #' @examples
# Purist=function(newexp){

#   # if(nrow(newexp)!= length(geneSymbols)){
#   #   stop("geneSymbols should be a vector of gene symbols exactly corresponding to each row of the newexp dataset")
#   # }
#   # expg=qutils::getUniqueGeneMat(newexp,geneSymbols,rowSds(as.matrix(newexp)))

#   .checkGeneExp(newexp)
#   inter=-6.815


#   puristcoef=data.frame(
#     a=c("GPR87","KRT6A","BCAR3","PTGES","ITGA3","C16orf74","S100A2","KRT5"),
#     b=c("REG4","ANXA10","GATA6","CLDN18","LGALS4","DDC","SLC40A1","CLRN3"),
#     w=c(1.994,2.031,1.618,0.922,1.059,0.929,2.505,0.485),stringsAsFactors = F)
#   puristg=c(puristcoef$a,puristcoef$b)


#   if(!all(puristg %in% rownames(newexp))){
#     stop(paste("Missing some genes:",paste(puristg[which(!puristg %in% rownames(newexp))],collapse=", ")))
#   }



#   prednum=setNames(rowSums(apply(puristcoef,1,function(p){(newexp[p[1],]>newexp[p[2],])*as.numeric(p[3])}))+inter,colnames(newexp))
#   preds=setNames(factor(c("classic","basal")[ (prednum>0)+1]),colnames(newexp))


#   data.frame(NumPurist=prednum,
#              Purist=preds,
#              row.names=colnames(newexp),stringsAsFactors = F)
# }
