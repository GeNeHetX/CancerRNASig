.getugm=function (m, g, w) 
{
    if (!all.equal(nrow(m), length(g), length(w))) {
        stop("nrow of m should be equal to lenght of g and w")
    }
    i = order(w, decreasing = T)
    oki = i[which(!duplicated(g[i]) & !g[i] %in% c("---", " ", 
        "", NA))]
    okm = m[oki, ]
    rownames(okm) = g[oki]
    okm
}


#' Title genesymUniqExp
#'
#' This function aims to obtain a gene expression matrix with a unique gene expression value
#' for each provided gene symbol. By default, the gene
#' 
#' @param newexp gene expression matrix/dataframe, sample in columns, gene in rows
#' @param geneSymbols gene symbols, a vector of same length as the number of rows in newex
#' @param scoreFunc a function taking a matrix and returning a score per line. this score is then
#'  used to select the highest scoring row if several correspond to the same gene symbol.
#'  With default rowSds, this means selecting the row with the highest variance.
#'
#' @return a data frame with row names as colnames of newexp, first column is purist subtype, second is purist score
#' @export
#'
#' @examples
genesymUniqExp = function(newexp,geneSymbols,scoreFunc=rowSds){
   .checkGeneExp(newexp)
    if(nrow(E)!=length(geneSymbols)){
        stop("Rows of newexp (expected to be genes) should have the same length as the proposed geneSymbols")
    }
    .getugm(E,geneSymbols,scoreFunc(E))

}

.checkGeneExp=function(newexp){


    if(!is.matrix(newexp)&!is.data.frame(newexp)){
        stop("newexp should be a matrix or dataframe")
    }
    E=as.matrix(newexp)
    if(ncol(E)<3){stop("newexp should have at least 3 columns")}
    if(nrow(E) < ncol(E)){
        warning("There are less rows than columns, genes are expected to be in lines (in the old bioinfomartics fashion)")
    }
}