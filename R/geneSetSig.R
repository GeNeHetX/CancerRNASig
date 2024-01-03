




#' Title coexpScore
#'
#' @param newexp A new expression dataset, should be gene symbol in rows
#' @param genesetlist Gene sets (signature or pathways) scores
#'
#' @return Gene set score
#' @export
#'
#' @examples
coexpScore=function(newexp,genesetlist){

    E=CancerRNASig:::.checkGeneExp(newexp)
    E=E[which(matrixStats::rowSds(E)> (10^-12)),]
    E <- t(base::scale(t(E), center=TRUE, scale = FALSE))

    CancerRNASig:::.checkCoregArgs(E, genesetlist)

    pp <- CancerRNASig:::.gesecaPreparePathways(E, genesetlist, minSize=1, 
    maxSize=(nrow(E)-1))
    genesetFiltered <- pp$filtered
    genesetSizes <- pp$sizes
    genesetNames <- names(genesetFiltered)

    genesetScores <- sapply(genesetFiltered, function(pa,E){
        return((matrixStats::colSums2(E[pa, , drop=FALSE])))
    }, E = E)

    genesetScores

}

.checkCoregArgs=function (E, genesetlist) 
{
    if (!is.list(genesetlist)) {
        stop("genesetlist should be a list with each element containing names of the stats argument")
    }
    if (is.null(rownames(E))) {
        stop("E rows should be named")
    }
    if (any(!is.finite(E))) {
        stop("Not all E values are finite numbers")
    }
    if (any(duplicated(rownames(E)))) {
        warning("There are duplicate gene names, geseca may produce unexpected results.")
    }
}


.gesecaPreparePathways <- function(E, pathways, minSize, maxSize){
    minSize <- max(minSize, 1)
    maxSize <- min(nrow(E) - 1, maxSize)

    pathwaysFiltered <- lapply(pathways, function(p) {unique(na.omit(match(p, rownames(E))))})
    pathwaysSizes <- sapply(pathwaysFiltered, length)

    toKeep <- which(minSize <= pathwaysSizes & pathwaysSizes <= maxSize)
    pathwaysFiltered <- pathwaysFiltered[toKeep]
    pathwaysSizes <- pathwaysSizes[toKeep]

    return(list(filtered=pathwaysFiltered,
                sizes=pathwaysSizes))
}


.evalOtherSig=function(g,gset,otherset){
    thisigscore=coregScore(E,list(a=setdiff(gset,c(otherset,g))))
    othesigscore=coregScore(E,list(b=otherset))
    blm=lm(E[g,]~othesigscore)
    alm=lm(E[g,]~thisigscore)

    balm=lm(E[g,]~othesigscore+thisigscore)
    residbalm=lm(blm$residuals~thisigscore)
    residablm=lm(alm$residuals~othesigscore)

    df=data.frame(gene=g,
        otheR2=summary(blm)$r.squared,
        othePv=summary(blm)$coefficients[2,4],
        gsetR2=summary(alm)$r.squared,
        bothR2=summary(balm)$r.squared,
        residotR2=summary(residbalm)$r.squared,
        residtoR2=summary(residablm)$r.squared,
        residPv=summary(residbalm)$coefficients[2,4],
        anova=anova(blm,balm)[2,6]
    )
    df$deltaR=df$gsetR2 -df$otheR2

    df
}

.twowaycomp=function(siga,sigb){
    siganob=setdiff(siga,sigb)
    sigbnoa=setdiff(sigb,siga)
    rbind( do.call(rbind,lapply(siganob,.evalOtherSig,gset=siganob,otherset=sigbnoa)),
   do.call(rbind,lapply(sigbnoa,.evalOtherSig,gset=sigbnoa,otherset=siganob)))
}



if(F){
    data("exampleExpressionMatrix")
    E=exampleExpressionMatrix

    cc=cor(t(E))
    hc=hclust(as.dist(1-cc))
    cl=cutree(hc,h=0.5)
    cli=tail(sort(table(cl)))[1:3]
    comclg=names(cl)[which(cl==names(cli)[1])]
    comclg2=names(cl)[which(cl==names(cli)[3])]
    cor(t(E[comclg,]))
    cor(t(E[comclg,]),t(E[comclg2,]))
    siga=sample(comclg,10)
    sigb=sample(comclg,10)
    
    same10=.twowaycomp(sample(comclg,10),sample(comclg,10))
    diff10=.twowaycomp(sample(comclg,10),sample(comclg2,10))
    boxplot(same10[,c("residotR2","residtoR2")])
    boxplot(diff10[,c("residotR2","residtoR2")])

    boxplot(list(
    same10=.twowaycomp(sample(comclg,10),sample(comclg,10))$residotR2,
    same10=.twowaycomp(sample(comclg,100),sample(comclg,100))$residotR2,
    diff10=.twowaycomp(sample(comclg,10),sample(comclg2,10))$residotR2,
    diff50=.twowaycomp(sample(comclg,50),sample(comclg2,50))$residotR2,
    diff100=.twowaycomp(sample(comclg,100),sample(comclg2,100))$residotR2
    ))


    
    boxplot(list(
    same10=.twowaycomp(sample(comclg,10),sample(comclg,10))$residtoR2,
    same10=.twowaycomp(sample(comclg,100),sample(comclg,100))$residtoR2,
    diff10=.twowaycomp(sample(comclg,10),sample(comclg2,10))$residtoR2,
    diff50=.twowaycomp(sample(comclg,50),sample(comclg2,50))$residtoR2,
    diff100=.twowaycomp(sample(comclg,100),sample(comclg2,100))$residtoR2
    ))


    boxplot(list(
    same10=.twowaycomp(sample(comclg,10),sample(comclg,10))$deltaR,
    same100=.twowaycomp(sample(comclg,100),sample(comclg,100))$deltaR,
    diff10=.twowaycomp(sample(comclg,10),sample(comclg2,10))$deltaR,
    diff50=.twowaycomp(sample(comclg,50),sample(comclg2,50))$deltaR,
    diff100=.twowaycomp(sample(comclg,100),sample(comclg2,100))$deltaR
    ))

       boxplot(list(
    same10=-log10(.twowaycomp(sample(comclg,10),sample(comclg,10))$anova),
    same10=-log10(.twowaycomp(sample(comclg,100),sample(comclg,100))$anova),
    diff10=-log10(.twowaycomp(sample(comclg,10),sample(comclg2,10))$anova),
    diff50=-log10(.twowaycomp(sample(comclg,50),sample(comclg2,50))$anova),
    diff100=-log10(.twowaycomp(sample(comclg,100),sample(comclg2,100))$anova)
    ))
    
    
    # g=siganob[1];  gset=siganob;    otherset=sigbnoa

   


    # system.time(mclapply(1:100,function(x) order(rnorm(n=1e6)),mc.cores=4))
    # system.time(BiocParallel::bplapply(1:100 , function(x) order(rnorm(n=1e6)), BPPARAM = MulticoreParam(workers = 4)))
 
}