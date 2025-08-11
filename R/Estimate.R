#' Title qEstimate
#'
#'
#' @return 
#' @export
#'
#' @examples
.estimate=function (ds, platform = c("illumina","affymetrix", "agilent"))
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
    gs <- as.matrix(CancerRNASig::estimategenes[, -1], dimnames = NULL)
    N.gs <- 2
    gs.names <- row.names(CancerRNASig::estimategenes)
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
    score.data <- rbind(score.data, estimate.score)
    rownames(score.data) <- c("StromalScore", "ImmuneScore","ESTIMATEScore")

    return(score.data)
}
