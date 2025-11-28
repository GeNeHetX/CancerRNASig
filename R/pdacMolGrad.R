#' Title projectMolGrad
#'
#' @description Project a transcriptomic dataset on the Pancreatic Adenocarcinoma Molecular Gradient
#' Will throw an error if there is only one sample (one column)
#'
#' @param newexp gene expression matrix or dataframe with gene in row(names) and samples in columns(names).
#' @param geneSymbols vector of gene symbols for the newexp dataset (simply set to rownames if the newexp is already in single values per gene symbols)
#' #' @param normalize Normalization (i.e. calibration) of the molecular gradient systems
#'
#' @return data frame of four projections based on the molecular gradients computed from four different types of expression datasets
#'
#' @details
#'
#'
#' @keywords internal
#'
#'
#' @examples
#' g <- rownames(CancerRNASig::molGradSys$PDX$gw)
#' projectMolGrad(matrix(rnorm(length(g) * 10), ncol = 10), g)
.projectMolGrad <- function(newexp, geneSymbols, normalize = c("newRef", "sameAsRef", "raw")) {
  data(molGradsys)
  normalize <- match.arg(normalize)
  if (nrow(newexp) != length(geneSymbols)) {
    stop("geneSymbols should be a vector of gene symbols exactly corresponding to each row of the newexp dataset")
  }
  # expg=qutils::getUniqueGeneMat(newexp,geneSymbols,matrixStats::rowSds(as.matrix(newexp)))
  projsL <- lapply(molGradsys, function(mg) {
    proj <- .internalProjection(newexp, mg, nlim = 4000, center = T, scale = T)
    if (is.na(proj[1])) {
      return(rep(NA, ncol(newexp)))
    }
    switch(normalize,
      newRef = {
        fproj <- ((proj - mean(proj)) / (3 * sd(proj)))
      },
      sameAsRef = {
        fproj <- ((proj - mg$avg) / (mg$sd))
      },
      raw = {
        fproj <- proj
      },
      {
        fproj <- proj
      }
    )
    return(fproj)
  })
  if (all(sapply(lapply(projsL, is.na), all))) {
    stop("geneSymbols should be a vector of Hugo gene symbols")
  }
  data.frame(do.call(cbind, projsL))
}

.internalProjection <- function(expg, sys, nlim = 4000, center = T, scale = T) {
  comg <- intersect(rownames(sys$gw), rownames(expg))
  if (length(comg) < nlim) {
    return(NA)
  } else {
    invs <- MASS::ginv(as.matrix(sys$gw[comg, ]))
    scexp <- scale(expg[comg, ], center = center, scale = scale)
    return((t(scexp) %*% t(invs))[, sys$k] * sys$dir)
  }
}
