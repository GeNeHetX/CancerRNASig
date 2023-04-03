#' Title qEstimate, default using Puleo et al components
#'
#' @param newexp A new expression dataset
#' @param ICAgw ICA gene weights (S matrix)
#' @param geneNormType Normalization of gene expression matrix (sc: sample scale is default)
#' @param projNormType Normalization of components (raw:  is default)
#' @param ming Minimum number of overlapping genes
#'
#' @return projected components
#' @export
#'
#' @examples
qEstimate=function (ds, platform = c("illumina","affymetrix", "agilent"
    ))
{
  require(estimate)

    platform <- match.arg(platform)
    # ds <- read.delim(input.ds, header = TRUE, sep = "\t", skip = 2,
    #     row.names = 1, blank.lines.skip = TRUE, as.is = TRUE,
    #     na.strings = "")
    # descs <- ds[, 1]
    # ds <- ds[-1]
    descs=NA
    row.names <- row.names(ds)
    names <- names(ds)
    dataset <- list(ds = ds, row.names = row.names, descs = descs,
        names = names)
    m <- data.matrix(dataset$ds)
    gene.names <- dataset$row.names
    sample.names <- dataset$names
    Ns <- length(m[1, ])
    Ng <- length(m[, 1])
    # temp <- strsplit(input.ds, split = "/")
    # s <- length(temp[[1]])
    # input.file.name <- temp[[1]][s]
    # temp <- strsplit(input.file.name, split = ".gct")
    # input.file.prefix <- temp[[1]][1]
    for (j in 1:Ns) {
        m[, j] <- rank(m[, j], ties.method = "average")
    }
    m <- 10000 * m/Ng
    gs <- as.matrix(estimate::SI_geneset[, -1], dimnames = NULL)
    N.gs <- 2
    gs.names <- row.names(estimate::SI_geneset)
    score.matrix <- matrix(0, nrow = N.gs, ncol = Ns)
    for (gs.i in 1:N.gs) {
        gene.set <- gs[gs.i, ]
        gene.overlap <- intersect(gene.set, gene.names)
        print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=",
            length(gene.overlap)))
        if (length(gene.overlap) == 0) {
            score.matrix[gs.i, ] <- rep(NA, Ns)
            next
        }
        else {
            ES.vector <- vector(length = Ns)
            for (S.index in 1:Ns) {
                gene.list <- order(m[, S.index], decreasing = TRUE)
                gene.set2 <- match(gene.overlap, gene.names)
                correl.vector <- m[gene.list, S.index]
                TAG <- sign(match(gene.list, gene.set2, nomatch = 0))
                no.TAG <- 1 - TAG
                N <- length(gene.list)
                Nh <- length(gene.set2)
                Nm <- N - Nh
                correl.vector <- abs(correl.vector)^0.25
                sum.correl <- sum(correl.vector[TAG == 1])
                P0 <- no.TAG/Nm
                F0 <- cumsum(P0)
                Pn <- TAG * correl.vector/sum.correl
                Fn <- cumsum(Pn)
                RES <- Fn - F0
                max.ES <- max(RES)
                min.ES <- min(RES)
                if (max.ES > -min.ES) {
                  arg.ES <- which.max(RES)
                }
                else {
                  arg.ES <- which.min(RES)
                }
                ES <- sum(RES)
                EnrichmentScore <- list(ES = ES, arg.ES = arg.ES,
                  RES = RES, indicator = TAG)
                ES.vector[S.index] <- EnrichmentScore$ES
            }
            score.matrix[gs.i, ] <- ES.vector
        }
    }
    score.data <- data.frame(score.matrix)
    names(score.data) <- sample.names
    row.names(score.data) <- gs.names
    estimate.score <- apply(score.data, 2, sum)
    if (platform != "affymetrix") {
        score.data <- rbind(score.data, estimate.score)
        rownames(score.data) <- c("StromalScore", "ImmuneScore",
            "ESTIMATEScore")
    }
    else {
        convert_row_estimate_score_to_tumor_purity <- function(x) {
            stopifnot(is.numeric(x))
            cos(0.6049872018 + 0.0001467884 * x)
        }
        est.new <- NULL
        for (i in 1:length(estimate.score)) {
            est_i <- convert_row_estimate_score_to_tumor_purity(estimate.score[i])
            est.new <- rbind(est.new, est_i)
            if (est_i >= 0) {
                next
            }
            else {
                message(paste(sample.names[i], ": out of bounds",
                  sep = ""))
            }
        }
        colnames(est.new) <- c("TumorPurity")
        estimate.t1 <- cbind(estimate.score, est.new)
        x.bad.tumor.purities <- estimate.t1[, "TumorPurity"] <
            0
        estimate.t1[x.bad.tumor.purities, "TumorPurity"] <- NA
        score.data <- rbind(score.data, t(estimate.t1))
        rownames(score.data) <- c("StromalScore", "ImmuneScore",
            "ESTIMATEScore", "TumorPurity")
    }
    return(score.data)
}
# qEstimate(tcgaexp)