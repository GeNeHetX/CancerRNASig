#' Title GemPred
#'
#' @description WARNING: This tool is intended for research purposes only. It is not intended for diagnostic use or clinical decision-making.
#' 
#' @param newexp gene expression matrix/dataframe, sample in columns, gene in rows
#' @param geneSymbols gene symbols, a vector of same length as the number of rows in newex
#'
#' @return a data frame with row names as colnames of newexp, first column is gempred score, second is gempred sensitivity conclusion
#'
#' @examples
#' # Example of how to call GemPred
#' # result <- GemPred(newexp = gene_expression_matrix, geneSymbols = gene_symbols)
#'
#' @export

GemPred = function(newexp, geneSymbols, q=0.75){
  
  newexp_GS = qutils::getUniqueGeneMat(newexp, geneSymbols, rowMeans(newexp))

  res_list <- lapply(colnames(newexp_GS), .lunch_GremPred, newexp_GS)
  names(res_list) <- colnames(newexp_GS)

  gempred_df <- do.call(rbind,res_list)

  # Compute sensitivity
  gempred_df$sensitivityGem <- "gp-"
  q3 <- quantile(gempred_df$gp2sc,probs = q)
  gempred_df[gempred_df$gp2sc>=q3,"sensitivityGem"] <- "gp+"
  gempred_df$sensitivityGem <- as.factor(gempred_df$sensitivityGem)
  res = gempred_df[,c("gp2sc","sensitivityGem")]
  
  return(res)
}

.lunch_GremPred<- function(name,counts){
    tmp_table <- counts[,c(name)]
    names(tmp_table) <- row.names(counts)
    pred = GemPred2(as.matrix(tmp_table))
    return(pred)
}

.getUGM=function (m, g, w)
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



GemPred1=function(dat){
  data(GWM)
  data(refAvg)
  w=c(-2.428483,0,0,1.716179,0,0)


  # expg=getUniqueGeneMat(dat$exp,dat$probeannot[rownames(dat$exp),dat$genecol],rowSds(as.matrix(dat$exp)))
  comg = intersect(rownames(dat),rownames(GWM))
  if(length(comg)<5000 ){
    return("Too few gene in common. (Gene Symbols are used)")
  }
  invs = MASS::ginv(as.matrix(GWM[comg, ]))
  v=as.matrix(dat[comg,])

  sscP=((t( invs %*% scale(v)) %*% w)[,1])*100
  rawP=(t( invs %*% (v)) %*% w)[,1]



  wwcomg=intersect(comg,rownames(refAvg))
  wwinvs = MASS::ginv(as.matrix(GWM[wwcomg, ]))

  vn=as.matrix(((dat[wwcomg,]-refAvg[wwcomg,"avg"])/refAvg[wwcomg,"sd"]))
  # vn=(((pancg[wwcomg,]-pancgRefAvg[wwcomg,"avg"])/pancgRefAvg[wwcomg,"sd"]))
  refedP=((t( wwinvs %*% ((as.matrix(dat[wwcomg,]-refAvg[wwcomg,"avg"])/refAvg[wwcomg,"sd"]))) %*% w)[,1])- -0.0786

  # c(rawP=rawP,
  data.frame(sscP=sscP, refedP=refedP,
             prob=dnorm(abs(apply(vn,2,mean))),
             setNames(data.frame(t((invs %*% v)[which(w !=0),])),c("senscomp","prolifcomp"))
  )

}

GemPred2=function(dat,useGeneSym=T,mingenes=500,compScale=F,preScale=F){
  data(GP2model)
  if(useGeneSym){
    gp2GW=CancerRNASig:::.getUGM(GP2model$ICA$S,GP2model$genes,apply(GP2model$ICA$S,1,max))
  }else{
    gp2GW=GP2model$ICA$S
  }
  comg2 = intersect(rownames(dat),rownames(gp2GW))

  if(length(comg2)<mingenes ){
    return(paste("Too few gene in common. (",length(comg2),")"))
  }

  invs2 = MASS::ginv(as.matrix(gp2GW[comg2, ]))
  v2=dat[comg2,]
  sscP2=t( invs2 %*% scale(v2,center=preScale,scale=preScale))
  sscP2=scale(sscP2,scale=compScale,center=compScale)
  colnames(sscP2)=colnames(gp2GW)

  gp2sc=unname(predict(GP2model$LM,data.frame(sscP2)))

  gp2Sr=crossprod(GP2model$simplified[comg2] , v2)[1,]

  res <- data.frame(gp2sc=gp2sc,
             refedP=factor(c("gp-","gp+")[1+(gp2sc > GP2model$cutoff)]),
             gp2singlegw=gp2Sr,
             PROJ=sscP2,
             row.names=colnames(dat)
  )

  res

}

GemPred_Lucie=function(dat,useGeneSym=T,mingenes=500){
  #data(GP2model_simple)
  if(useGeneSym){
    gp2GW=GP2model_simple$GeneSymbol
  }else{
    gp2GW=GP2model_simple$ENG_ID
  }
  comg2 = intersect(rownames(dat),rownames(gp2GW))

  if(length(comg2)<mingenes ){
    return(paste("Too few gene in common. (",length(comg2),")"))
  }
  v2=dat[comg2,]
  gp2Sr=crossprod(GP2model_simple$simplified[comg2] , v2)[1,]

  res <- data.frame(gp2sc=gp2Sr,
             row.names=colnames(dat)
  )

  res

}
