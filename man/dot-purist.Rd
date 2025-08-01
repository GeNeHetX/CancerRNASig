% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Purist.R
\name{.purist}
\alias{.purist}
\title{Title purist}
\usage{
.purist(newexp, geneSymbols)
}
\arguments{
\item{newexp}{gene expression matrix/dataframe, sample in columns, gene in rows}

\item{geneSymbols}{gene symbols, a vector of same length as the number of rows in newex}
}
\value{
a data frame with row names as colnames of newexp, first column is purist subtype, second is purist score
}
\description{
Title purist
}
\keyword{internal}
