
#' Title Simple ICA projection, default using Puleo et al components
#'
#' @param newexp A new expression dataset, should be gene symbol in rows
#' @param ICAgw ICA gene weights (S matrix)
#' @param geneNormType Normalization of gene expression matrix (sc: sample scale is default)
#' @param projNormType Normalization of components (raw:  is default)
#' @param ming Minimum number of overlapping genes
#'
#' @return projected components
#' @export
#'
#' @examples
qProjICA=function(newexp,ICAgw=CancerRNASig:::puleoICAgw,geneNormType="sc",projNormType="raw",ming=500){

  comg = intersect(rownames(newexp), rownames(ICAgw))

  if(length(comg)<ming){
    stop("Too few genes in common in the expression dataset and the ICA weights")
  }

  scexp = CancerRNASig:::.qNormalize(newexp[comg, ],type=geneNormType)

  invs = MASS::ginv(as.matrix(ICAgw[comg,]))

  proj=t(CancerRNASig:::.qNormalize(invs %*%  scexp,projNormType))

  colnames(proj)=colnames(ICAgw)
  proj

}


.qNumerote=function(n,l=max(nchar(1:n))) {
  str_pad(1:n,width=l,side="left",pad=0)
}

.qNormalize=function(x,type){
  switch(type,
         raw={x},
         gsc={t(scale(t(x)))},
         gc={t(scale(t(x),scale=F))},
         sc={(scale((x),scale=F))},
         ssc={(scale((x)))}
         ,x
  )
}



.qICA=function(X,k,maxiter = 10^6,eps = 10^-6){
  resICA=NULL
  try({
    resICA <- JADE::JADE(X, n.comp = k, maxiter = maxiter, eps = eps)
    rownames(resICA$A) <- colnames(X)
    rownames(resICA$S) <- rownames(X)
    colnames(resICA$S)=paste("ICA",CancerRNASig:::.qNumerote(ncol(resICA$S)),sep="")
    colnames(resICA$A)=paste("ICA",CancerRNASig:::.qNumerote(ncol(resICA$A)),sep="")
  })
  resICA
}






.qProjICA.ds=function(icarez,dataset,geneNormType="sc",projNormType="raw"){

  # expg=getUniqueGeneMat(dataset$exp,dataset$probeannot[rownames(dataset$exp),dataset$genecol],rowSds(as.matrix(dataset$exp)))
   comg = intersect(rownames(icarez$S), rownames(dataset$exp))
  scexp = CancerRNASig:::.qNormalize(dataset$exp[comg, ],type=geneNormType)


  invs = MASS::ginv(as.matrix(icarez$S[comg,]))

  proj=t(CancerRNASig:::.qNormalize(invs %*%  scexp,projNormType))

  # proj = scale(t(t((t(scexp) %*% t(invs)))))

  colnames(proj)=colnames(icarez$S)
  proj

}

