% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GemPred.R
\name{.gemPred}
\alias{.gemPred}
\title{Title GemPred}
\usage{
.gemPred(newexp, geneSymbols, q = 0.75)
}
\arguments{
\item{newexp}{gene expression matrix/dataframe, sample in columns, gene in rows}

\item{geneSymbols}{gene symbols, a vector of same length as the number of rows in newex}
}
\value{
a data frame with row names as colnames of newexp, first column is gempred score, second is gempred sensitivity conclusion
}
\description{
WARNING: This tool is intended for research purposes only. It is not intended for diagnostic use or clinical decision-making.
}
\examples{
# Example of how to call GemPred
# result <- GemPred(newexp = gene_expression_matrix, geneSymbols = gene_symbols)

}
\keyword{internal}
