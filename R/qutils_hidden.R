# Hidden compatibility helpers ported from qutils.
# These are intentionally not exported from the package namespace.

qNumerote <- function(n, l = max(nchar(1:n))) {
  stringr::str_pad(1:n, width = l, side = "left", pad = 0)
}

.qNumerote <- qNumerote

qnarm <- narm <- function(x) x[which(!is.na(x))]

qemptyrm <- emptyrm <- function(x) x[which(x != "")]

quantNorm2ref <- function(x, ref) {
  if (nrow(ref) != nrow(x)) {
    stop("both datsets need to have same number of rows")
  }

  xr <- apply(x, 2, rank, ties.method = "min")
  refs <- data.frame(apply(ref, 2, sort))
  refm <- rowMeans(as.matrix(refs))

  index_to_mean <- function(my_index, my_mean) {
    my_mean[my_index]
  }

  df_final <- apply(xr, 2, index_to_mean, my_mean = refm)
  dimnames(df_final) <- dimnames(x)
  df_final
}

getUniqueGeneMat <- function(m, g, w) {
  if (!all.equal(nrow(m), length(g), length(w))) {
    stop("nrow of m should be equal to lenght of g and w")
  }

  i <- order(w, decreasing = TRUE)
  oki <- i[which(!duplicated(g[i]) & !g[i] %in% c("---", " ", "", NA))]

  okm <- m[oki, ]
  rownames(okm) <- g[oki]
  okm
}

getUniqueGeneVec <- function(m, g, w) {
  if (!all.equal(length(m), length(g), length(w))) {
    stop("length of m should be equal to lenght of g and w")
  }

  i <- order(w, decreasing = TRUE)
  oki <- i[which(!duplicated(g[i]) & !g[i] %in% c("---", " ", "", NA))]

  okm <- m[oki]
  names(okm) <- g[oki]
  okm
}

autopathrenam <- function(path) {
  gsub("_Homo sapiens_.+", "", gsub("^[^~]+~", "", path))
}

qfgsea <- function(signedcor, pathways, n = 1000, thresh = 0.01) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required for qfgsea().")
  }
  set.seed(1)
  xgsea <- fgsea::fgseaSimple(pathways, signedcor, nperm = n, nproc = 1)
  xgsea <- xgsea[which(xgsea$pval < thresh), ]
  xgsea[order(-abs(xgsea$NES)), ]
}

qGseaTable <- function(pathways, stats, fgseaRes, gseaParam = 1,
                       colwidths = c(5, 3, 0.8, 1.2, 1.2), rename = NULL) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required for qGseaTable().")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for qGseaTable().")
  }
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package 'gridExtra' is required for qGseaTable().")
  }

  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  n <- length(statsAdj)

  pathways <- lapply(pathways, function(p) {
    unname(as.vector(na.omit(match(p, names(statsAdj)))))
  })

  if (!is.null(rename) && length(rename) == length(pathways)) {
    names(rename) <- names(pathways)
  } else {
    rename <- names(pathways)
    names(rename) <- names(pathways)
  }

  ps <- lapply(names(pathways), function(pn) {
    p <- pathways[[pn]]
    annotation <- fgseaRes[match(pn, fgseaRes$pathway), ]
    list(
      grid::textGrob(rename[pn], just = "right", x = grid::unit(0.95, "npc")),
      ggplot2::ggplot() +
        ggplot2::geom_segment(
          ggplot2::aes(x = p, xend = p, y = 0, yend = statsAdj[p], col = p),
          linewidth = 0.2
        ) +
        ggplot2::scale_x_continuous(limits = c(0, length(statsAdj)), expand = c(0, 0)) +
        ggplot2::scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
        ggplot2::xlab(NULL) + ggplot2::ylab(NULL) +
        ggplot2::scale_colour_gradient2(
          midpoint = n / 2, low = "#ca0020", mid = "#f7f7f7", high = "#0571b0"
        ) +
        ggplot2::theme(
          panel.background = ggplot2::element_blank(),
          legend.position = "none",
          axis.line = ggplot2::element_blank(),
          axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank(),
          axis.title = ggplot2::element_blank(),
          plot.margin = grid::unit(rep(0, 4), "null"),
          panel.spacing = grid::unit(rep(0, 4), "null")
        ),
      grid::textGrob(sprintf("%.2f", annotation$NES)),
      grid::textGrob(sprintf("%.1e", annotation$pval)),
      grid::textGrob(sprintf("%.1e", annotation$padj))
    )
  })

  rankPlot <- ggplot2::ggplot() +
    ggplot2::geom_blank() +
    ggplot2::scale_x_continuous(limits = c(0, length(statsAdj)), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
    ggplot2::xlab(NULL) + ggplot2::ylab(NULL) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      plot.margin = grid::unit(c(0, 0, 0.5, 0), "npc"),
      panel.spacing = grid::unit(c(0, 0, 0, 0), "npc")
    )

  grobs <- c(
    lapply(c("Pathway", "Gene ranks", "NES", "pval", "padj"), grid::textGrob),
    unlist(ps, recursive = FALSE),
    list(grid::nullGrob(), rankPlot, grid::nullGrob(), grid::nullGrob(), grid::nullGrob())
  )
  grobsToDraw <- rep(colwidths != 0, length(grobs) / length(colwidths))
  gridExtra::grid.arrange(
    grobs = grobs[grobsToDraw],
    ncol = sum(colwidths != 0),
    widths = colwidths[colwidths != 0]
  )
}

qGseaPlot <- function(pathwayName, stats, fgseaRez, pathwayDatabase, main = NULL,
                      signifN = 3, gseaParam = 1, doprint = TRUE) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required for qGseaPlot().")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for qGseaPlot().")
  }

  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathwayDatabase[[pathwayName]], names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- fgsea:::calcGseaStat(statsAdj, selectedStats = pathway, returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms)) / 8

  if (abs(max(tops)) > abs(min(bottoms))) {
    upper <- max(tops) + 0.15
    lower <- 0 - diff / 2
    laby <- upper - 0.05
  } else {
    lower <- min(bottoms) - 0.1
    upper <- diff / 2 + 0.1
    laby <- upper - 0.05
  }

  fgsea_row <- fgseaRez[match(pathwayName, fgseaRez$pathway), , drop = FALSE]

  g <- ggplot2::ggplot(toPlot, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") +
    ggplot2::geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, colour = "black") +
    ggplot2::geom_line(color = "green") +
    ggplot2::theme_bw() +
    ggplot2::geom_segment(
      data = data.frame(x = pathway),
      ggplot2::aes(x = x, y = -diff / 2, xend = x, yend = diff / 2),
      linewidth = 0.2
    ) +
    ggplot2::geom_segment(
      data = data.frame(x = pathway),
      ggplot2::aes(x = x, y = -diff / 2, xend = x, yend = diff / 2, colour = x),
      linewidth = 0.2
    ) +
    ggplot2::scale_colour_gradient2(
      midpoint = n / 2, low = "#ca0020", mid = "#f7f7f7", high = "#0571b0"
    ) +
    ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(x = "rank", y = "enrichment score") +
    ggplot2::guides(colour = "none") +
    ggplot2::ylim(c(lower, upper)) +
    ggplot2::xlim(0, n + 500) +
    ggplot2::annotate(
      "text", x = n / 2, y = laby,
      label = paste(
        "NES:", signif(fgsea_row$NES, signifN),
        "  p-value:", signif(fgsea_row$pval, signifN)
      )
    )

  if (!is.null(main)) {
    g <- g + ggplot2::ggtitle(main)
  }

  g <- g + ggplot2::theme_minimal()
  if (doprint) {
    print(g)
    invisible(g)
  } else {
    g
  }
}

zhmcol <- function(w = 1) circlize::colorRamp2((-1:1) * w, c("blue", "white", "red"))
ahmcol <- function(w = 1) circlize::colorRamp2((0:1) * w, c("white", "red"))
mhmcol <- function(w = 1) circlize::colorRamp2((0:1) * w, c("blue", "yellow"))

CCCP <- function(projectTitle,
                 data,
                 weights = NULL,
                 meth.Clustering = c("ward.D2", "complete", "average"),
                 meth.Distance = c("pearson", "euclidean", "spearman"),
                 nbootstrap = 1000,
                 pGenesProbes = 0.8,
                 pSamples = 0.8,
                 qselectd = exp(seq(log(0.01), log(0.5), length.out = 10)) * nrow(data),
                 Kmax = 10,
                 finalLinkage = "complete",
                 ANNOT = NULL) {
  if (!requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
    stop("Package 'ConsensusClusterPlus' is required for CCCP().")
  }
  if (!requireNamespace("NMF", quietly = TRUE)) {
    stop("Package 'NMF' is required for CCCP().")
  }
  if (!requireNamespace("SAGx", quietly = TRUE)) {
    stop("Package 'SAGx' is required for CCCP().")
  }
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Package 'cluster' is required for CCCP().")
  }

  qselectd <- ceiling(qselectd)
  dir.create(projectTitle, showWarnings = FALSE)
  oldwd <- getwd()
  setwd(projectTitle)
  on.exit(setwd(oldwd), add = TRUE)

  dat.norm <- as.matrix(data)
  if (is.null(weights)) {
    weights <- matrixStats::rowMads(dat.norm)
    names(weights) <- rownames(dat.norm)
  }

  grDevices::pdf("ALLCONSENSUS.pdf")
  on.exit(grDevices::dev.off(), add = TRUE)

  customCCCPrun <- parallel::mclapply(meth.Distance, function(dis) {
    arun <- parallel::mclapply(meth.Clustering, function(linkage) {
      arun <- lapply(qselectd, function(aq) {
        titre <- paste(projectTitle, aq, dis, linkage, sep = "_")
        if (aq > 0 && aq <= 1) {
          psok <- names(weights)[which(weights >= stats::quantile(weights, probs = aq))]
        }
        if (aq > 1 && aq <= nrow(dat.norm)) {
          psok <- rownames(dat.norm)[which(rank(weights) >= (length(weights) - (aq - 1)))]
        }
        if (aq == 0) {
          psok <- rownames(dat.norm)
        }

        ccp <- ConsensusClusterPlus::ConsensusClusterPlus(
          dat.norm[psok, ],
          maxK = Kmax, reps = nbootstrap, pItem = pSamples, pFeature = pGenesProbes,
          plot = NULL, clusterAlg = "hc", verbose = TRUE, title = titre,
          distance = dis, innerLinkage = linkage,
          finalLinkage = finalLinkage, corUse = "pairwise.complete.obs"
        )
        lapply(ccp[2:Kmax], function(x) x$consensusMatrix)
      })

      lapply(1:(Kmax - 1), function(k) {
        X <- Reduce("+", lapply(arun, function(x) x[[k]])) / length(qselectd)
        rownames(X) <- colnames(dat.norm)
        colnames(X) <- colnames(dat.norm)
        X
      })
    })

    lapply(1:(Kmax - 1), function(k) {
      X <- Reduce("+", lapply(arun, function(x) x[[k]])) / length(meth.Clustering)
      rownames(X) <- colnames(dat.norm)
      colnames(X) <- colnames(dat.norm)
      X
    })
  })

  cccpConsMat <- lapply(1:(Kmax - 1), function(k) {
    X <- Reduce("+", lapply(customCCCPrun, function(x) x[[k]])) / length(meth.Distance)
    rownames(X) <- colnames(dat.norm)
    colnames(X) <- colnames(dat.norm)
    X
  })

  grDevices::pdf("consensusClusterings.pdf")
  on.exit(grDevices::dev.off(), add = TRUE)
  lapply(seq_along(cccpConsMat), function(i) {
    if (is.null(ANNOT)) {
      graphics::heatmap(
        cccpConsMat[[i]][colnames(data), colnames(data)],
        scale = "none",
        col = grDevices::colorRampPalette(c("white", "blue"))(32),
        main = paste("Consensus clustering k=", i + 1)
      )
    } else {
      X <- data.frame(cccpConsMat[[i]], stringsAsFactors = FALSE)
      rownames(X) <- rownames(cccpConsMat[[i]])
      colnames(X) <- colnames(cccpConsMat[[i]])
      cit.heatmap(
        X, scale = "none", heatmapcolors = grDevices::colorRampPalette(c("white", "blue"))(32),
        colclust.hclust = cit.hclust(dx = as.dist(1 - cccpConsMat[[i]]), method = "complete"), colclust.k = i + 1,
        rowclust.hclust = cit.hclust(dx = as.dist(1 - cccpConsMat[[i]]), method = "complete"), rowclust.k = i + 1,
        colannot = ANNOT[colnames(cccpConsMat[[i]]), ],
        lw = c(12, 7, 1, 16), lh = c(3, 10, 1, 15),
        test = "chisq.test"
      )
    }
  })

  cccpConsHierarchy <- lapply(cccpConsMat, function(x) stats::hclust(stats::as.dist(1 - x), method = finalLinkage))
  cccpConsPart <- lapply(1:(Kmax - 1), function(k) stats::cutree(cccpConsHierarchy[[k]], k + 1))

  save(cccpConsPart, file = file.path("cccpConsensusPartition.RData"))
  save(cccpConsHierarchy, file = file.path("cccpConsensusHierarchy.RData"))
  save(cccpConsMat, file = file.path("cccpConsensusCoocurenceMatrix.RData"))

  q <- qselectd[length(qselectd)]
  if (q > 0 && q <= 1) qkeep <- names(weights)[which(weights >= stats::quantile(weights, probs = q))]
  if (q > 1 && q <= nrow(dat.norm)) qkeep <- rownames(dat.norm)[which(rank(weights) >= (length(weights) - (q - 1)))]
  if (q == 0) qkeep <- rownames(dat.norm)

  subdat.norm <- as.matrix(dat.norm[qkeep, ])

  cccpmetrics_subdata <- .cit.clustMetric(subdat.norm, cccpConsMat, cccpConsHierarchy, cccpConsPart, 2:Kmax, plot = TRUE, pdf.name = "CCCP.clustMetric_subdata.pdf")
  cccpmetrics_alldata <- .cit.clustMetric(dat.norm, cccpConsMat, cccpConsHierarchy, cccpConsPart, 2:Kmax, plot = TRUE, pdf.name = "CCCP.clustMetric_alldata.pdf")

  list(
    cccpConsPart = cccpConsPart,
    cccpConsHierarchy = cccpConsHierarchy,
    cccpConsMat = cccpConsMat,
    cccpmetrics_subdata = cccpmetrics_subdata,
    cccpmetrics_alldata = cccpmetrics_alldata
  )
}

.cit.clustMetric <- function(expdata, consensuscooccurencematrix, hierarchy, partitions,
                             ks = 2:ifelse(is.list(hierarchy), length(hierarchy) + 1, length(unique(partitions))),
                             plot = TRUE, pdf.name = "clustMetric.pdf", doGAPstat = FALSE) {
  if (!requireNamespace("NMF", quietly = TRUE)) {
    stop("Package 'NMF' is required for .cit.clustMetric().")
  }
  if (!requireNamespace("SAGx", quietly = TRUE)) {
    stop("Package 'SAGx' is required for .cit.clustMetric().")
  }
  originalDist <- stats::dist(t(expdata))
  if (plot) grDevices::pdf(pdf.name)

  on.exit(if (plot) grDevices::dev.off(), add = TRUE)

  if (is.list(consensuscooccurencematrix) &&
      length(consensuscooccurencematrix) == length(hierarchy) &&
      length(consensuscooccurencematrix) == length(partitions)) {
    rez <- data.frame(
      coph = sapply(hierarchy, function(h) cor(as.numeric(stats::cophenetic(h)), as.numeric(originalDist))),
      disp = sapply(consensuscooccurencematrix, NMF::dispersion),
      silhouettes = sapply(partitions, function(part) {
        p <- as.numeric(as.factor(part))
        names(p) <- names(part)
        sil <- cluster::silhouette(p, dist = originalDist)
        summary(sil)$avg.width
      })
    )
    if (doGAPstat) {
      rez$gap <- sapply(partitions, function(part) {
        try(SAGx::gap(t(expdata), part[colnames(expdata)])[1])
      })
    }
    if (plot) {
      if (doGAPstat) {
        plot(1:length(rez$gap), rez$gap, xlab = "Cluster number", ylab = "Gap statistics",
             main = "Gap Statistics", type = "b", axes = FALSE, pch = 16)
        axis(2, las = 2)
        axis(1, at = 1:length(rez$gap), labels = paste("K", 2:(length(rez$gap) + 1), sep = ""))
      }
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        y <- as.data.frame(as.table(as.matrix(rez)), stringsAsFactors = FALSE)
        colnames(y) <- c("K", "Metric", "value")
        if (!is.factor(y$K)) {
          y$K <- factor(paste("K=", y$K + 1, sep = ""), levels = paste("K=", ks, sep = ""))
        }
        print(ggplot2::ggplot(y, ggplot2::aes(x = K, y = value, col = Metric)) +
                ggplot2::geom_point() + ggplot2::theme_classic())
      } else {
        matplot(rez, type = "l", lwd = 3, x = 2:Kmax, xlab = "K", ylab = "Metric",
                ylim = c(min(as.matrix(rez)), 1.1))
        legend("topright", colnames(rez), col = 1:ncol(rez), horiz = TRUE, lty = 1:ncol(rez))
      }
      plot(1:length(rez$coph), rez$coph, xlab = "Cluster number", ylab = "Cophenetic correlation",
           main = "Cophenetic coefficient", type = "b", axes = FALSE, pch = 16)
      axis(2, las = 2)
      axis(1, at = 1:length(rez$coph), labels = paste("K", 2:(length(rez$coph) + 1), sep = ""))
      plot(1:length(rez$disp), rez$disp, xlab = "Cluster number", ylab = "Dispersion",
           main = "Dispersion", type = "b", axes = FALSE, pch = 16)
      axis(2, las = 2)
      axis(1, at = 1:length(rez$disp), labels = paste("K", 2:(length(rez$disp) + 1), sep = ""))
      plot(1:length(rez$silhouettes), rez$silhouettes, xlab = "Cluster number", ylab = "Mean silhouette",
           main = "Mean silhouettes", type = "b", axes = FALSE, pch = 16)
      axis(2, las = 2)
      axis(1, at = 1:length(rez$silhouettes), labels = paste("K", 2:(length(rez$silhouettes) + 1), sep = ""))
      N <- length(consensuscooccurencematrix)
      colv <- rainbow(N + 1)
      cdf <- .cit.computeCDF(consensuscooccurencematrix[[1]])
      plot(as.numeric(names(cdf)), cdf, xlab = "Consensus index value", ylab = "CDF",
           col = colv[1], type = "l", lwd = 2.5, main = "Cumulative Distribution Function")
      for (i in 2:N) {
        tmp <- .cit.computeCDF(consensuscooccurencematrix[[i]])
        lines(as.numeric(names(tmp)), tmp, col = colv[i], lwd = 2.5)
      }
      legend("topleft", legend = paste("K=", 2:(N + 1)), fill = colv[1:N], bty = "n")
      sapply(partitions, function(part) {
        p <- as.numeric(as.factor(part))
        names(p) <- names(part)
        sil <- cluster::silhouette(p, dist = originalDist)
        plot(sil)
        summary(sil)$avg.width
      })
    }
  } else {
    sil <- cluster::silhouette(partitions, dist = originalDist)
    rez <- data.frame(
      coph = cor(stats::cophenetic(hierarchy), originalDist),
      disp = NMF::dispersion(consensuscooccurencematrix),
      silhouettes = summary(sil)$avg.width
    )
    if (doGAPstat) rez$gap <- try(SAGx::gap(t(expdata), partitions[names(expdata)])[1])
    if (plot) {
      cdf <- .cit.computeCDF(consensuscooccurencematrix)
      plot(as.numeric(names(cdf)), cdf, xlab = "Consensus index value", ylab = "CDF",
           col = "red", type = "l", lwd = 2.5, main = "Cumulative Distribution Function")
      plot(sil)
    }
  }

  rownames(rez) <- paste("k=", ks, sep = "")
  rez
}

.cit.computeCDF <- function(coocurencematrix) {
  if (nrow(coocurencematrix) != ncol(coocurencematrix)) {
    stop("A co-classification (s x s) matrix is required")
  }
  consind <- seq(0, 1, 0.05)
  coocurencematrix <- as.matrix(coocurencematrix)
  upmat <- coocurencematrix[upper.tri(coocurencematrix)]
  N <- length(upmat)
  sapply(consind, function(ci) sum(upmat <= ci) / N)
}

rowMeanCorrectedSD <- function(X) {
  m <- rowMeans(X)
  s <- matrixStats::rowSds(X)
  r <- rank(matrixStats::rowMins(cbind(rank(m), rank(s)))) / nrow(X)
  names(r) <- rownames(X)
  r
}

qlimma <- function(expressiondata, samplesgrp1, samplesgrp2, thresh = 0.01, namesonly = TRUE, addWillcox = FALSE) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for qlimma().")
  }
  samplesgrp1 <- intersect(samplesgrp1, colnames(expressiondata))
  samplesgrp2 <- intersect(samplesgrp2, colnames(expressiondata))
  expressiondata <- expressiondata[, c(samplesgrp2, samplesgrp1)]
  classfactor <- factor(c(rep("a", length(samplesgrp2)), rep("b", length(samplesgrp1))))
  design <- model.matrix(~classfactor)
  colnames(design)[2] <- "avsb"
  contrast.matrix <- limma::makeContrasts("a-b", levels = classfactor)
  fit <- limma::lmFit(expressiondata, design)
  fit2 <- limma::eBayes(fit)
  fulltab <- limma::topTable(fit2, adjust = "BH", number = nrow(expressiondata))
  if (namesonly) {
    rownames(fulltab)[which(fulltab$adj.P.Val < thresh)]
  } else {
    X <- fulltab[which(fulltab$adj.P.Val < thresh), ]
    if (addWillcox) {
      X$Wilcoxon.Pvalue <- unlist(parallel::mclapply(rownames(X), function(prob) {
        stats::wilcox.test(
          as.numeric(expressiondata[prob, samplesgrp1]),
          as.numeric(expressiondata[prob, samplesgrp2])
        )$p.value
      }))
    }
    X
  }
}

fpDrawSimpleCI <- function(lower_limit, estimate, upper_limit, size = 1, y.offset = 0.5,
                           clr.line, clr.marker, lwd, lty = 1, vertices, vertices.height = 0.1,
                           ...) {
  if (is.na(lower_limit) || is.na(estimate) || is.na(upper_limit)) {
    return()
  }
  forestplot:::prFpDrawLine(
    lower_limit = lower_limit, upper_limit = upper_limit,
    clr.line = clr.line, lwd = lwd, lty = lty, y.offset = y.offset,
    vertices = vertices, vertices.height = vertices.height
  )
  box <- grid::convertX(grid::unit(estimate, "native"), "npc", valueOnly = TRUE)
  skipbox <- box < 0 || box > 1
  if (!skipbox) {
    h <- grid::unit(0.2, "snpc")
    w <- grid::unit(0.01, "snpc")
    grid::grid.rect(
      x = grid::unit(estimate, "native"), y = y.offset,
      width = w, height = h,
      gp = grid::gpar(fill = clr.marker, col = clr.marker)
    )
  }
}

.cipaste <- function(xv) {
  xv <- signif(xv, 3)
  paste0(xv[1], " (", xv[2], "-", xv[3], ")")
}

multivariateForest <- function(acoxph, signifnum = 3, ...) {
  if (!requireNamespace("forestplot", quietly = TRUE)) {
    stop("Package 'forestplot' is required for multivariateForest().")
  }
  sumcox <- summary(acoxph)
  tabletext <- list(
    c(NA, rownames(sumcox$conf.int)),
    c("HR (95% CI)", apply(sumcox$conf.int[, c(1, 3, 4)], 1, .cipaste)),
    c("P-value (Wald)", signif(sumcox$coefficients[, 5], signifnum))
  )
  toplotvals <- rbind(rep(NA, 3), sumcox$conf.int[, c(1, 3, 4)])
  forestplot::forestplot(
    tabletext, toplotvals, align = rep("l", 5), xlog = TRUE,
    is.summary = c(TRUE, rep(FALSE, nrow(sumcox$conf.int))),
    hrzl_lines = grid::gpar(col = "#444444"), graph.pos = 2,
    fn.ci_norm = fpDrawSimpleCI, ...
  )
}

univariateForest <- function(listCoxph, signifnum = 3, printArgHelp = FALSE, ...) {
  if (!requireNamespace("forestplot", quietly = TRUE)) {
    stop("Package 'forestplot' is required for univariateForest().")
  }
  if (printArgHelp) {
    print(paste(
      "clip: upper and lower limit of HR to plots (adds arrow if necessary)",
      "xlab for label of plot",
      "graph.pos : column to put plot of HR",
      "hrzl_lines = list('2'=gpar(col='#444444')) : to add a line under first row    ",
      "col = fpColors(box='royalblue',line='darkblue') : to color plots",
      "w", sep = "\n"
    ))
  }

  tabletext <- list(
    c(NA, names(listCoxph)),
    c("HR (95% CI)", sapply(listCoxph, function(x) {
      xx <- summary(x)
      paste0(signif(xx$conf.int[1, 1], signifnum), " (", signif(xx$conf.int[1, 3], signifnum), "-", signif(xx$conf.int[1, 4], signifnum), ")")
    })),
    c("P-value (Wald)", sapply(listCoxph, function(x) {
      xx <- summary(x)
      as.character(signif(xx$waldtest[3], signifnum))
    })),
    c("P-value (Score)", sapply(listCoxph, function(x) {
      xx <- summary(x)
      as.character(signif(xx$sctest[3], signifnum))
    }))
  )

  toplotvals <- rbind(rep(NA, 3), t(sapply(listCoxph, function(x) summary(x)$conf.int[1, c(1, 3, 4)])))
  forestplot::forestplot(
    tabletext, toplotvals, align = rep("l", 5), xlog = TRUE,
    is.summary = c(TRUE, rep(FALSE, length(listCoxph))),
    hrzl_lines = grid::gpar(col = "#444444"), graph.pos = 2,
    fn.ci_norm = fpDrawSimpleCI, ...
  )
}

dualunivariateForest <- function(listCoxph1, listCoxph2, signifnum = 3, cols = c("green", "brown"), ...) {
  if (!requireNamespace("forestplot", quietly = TRUE)) {
    stop("Package 'forestplot' is required for dualunivariateForest().")
  }
  if (any(names(listCoxph1) != names(listCoxph2))) {
    stop("not same names of list of cox")
  }

  tabletext <- list(
    c(NA, names(listCoxph1)),
    c("HR (95% CI)", sapply(names(listCoxph1), function(x) {
      xx1 <- summary(listCoxph1[[x]])
      xx2 <- summary(listCoxph2[[x]])
      paste0(
        signif(xx1$conf.int[1, 1], signifnum), " (", signif(xx1$conf.int[1, 3], signifnum), "-", signif(xx1$conf.int[1, 4], signifnum), ")\n",
        signif(xx2$conf.int[1, 1], signifnum), " (", signif(xx2$conf.int[1, 3], signifnum), "-", signif(xx2$conf.int[1, 4], signifnum), ")"
      )
    })),
    c("P-value (Wald)", sapply(names(listCoxph1), function(x) {
      xx1 <- summary(listCoxph1[[x]])
      xx2 <- summary(listCoxph2[[x]])
      paste0(as.character(signif(xx1$waldtest[3], signifnum)), "\n", as.character(signif(xx2$waldtest[3], signifnum)))
    })),
    c("P-value (Score)", sapply(names(listCoxph1), function(x) {
      xx1 <- summary(listCoxph1[[x]])
      xx2 <- summary(listCoxph2[[x]])
      paste0(as.character(signif(xx1$sctest[3], signifnum)), "\n", as.character(signif(xx2$sctest[3], signifnum)))
    }))
  )

  toplotvals1 <- rbind(rep(NA, 3), t(sapply(listCoxph1, function(x) summary(x)$conf.int[1, c(1, 3, 4)])))
  toplotvals2 <- rbind(rep(NA, 3), t(sapply(listCoxph2, function(x) summary(x)$conf.int[1, c(1, 3, 4)])))

  forestplot::forestplot(
    tabletext,
    mean = cbind(toplotvals1[, 1], toplotvals2[, 1]),
    lower = cbind(toplotvals1[, 2], toplotvals2[, 2]),
    upper = cbind(toplotvals1[, 3], toplotvals2[, 3]),
    xlog = TRUE, align = rep("l", 5),
    is.summary = c(TRUE, rep(FALSE, length(listCoxph1))),
    hrzl_lines = grid::gpar(col = "#444444"), graph.pos = 2,
    fn.ci_norm = fpDrawSimpleCI,
    col = forestplot::fpColors(box = cols), ...
  )
}
