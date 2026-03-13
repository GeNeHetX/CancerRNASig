#' Datasets included in CancerRNASig
#'
#' Collection of datasets used for RNA signature prediction and
#' tumor microenvironment characterization.
#'
#' @format Various objects used internally by the package.
#' @source Literature-derived gene signatures and models.
#' @name CancerRNASig_data
NULL

#' @rdname CancerRNASig_data
#' @docType data
"GP2model_simple"

#' @rdname CancerRNASig_data
#' @docType data
"molGradsys"

#' @rdname CancerRNASig_data
#' @docType data
"signatures"

#' @rdname CancerRNASig_data
#' @docType data
"t_GemPred_Biopsy"


#' Import required base R functions
#'
#' These imports make certain functions from stats and utils visible to R CMD check.
#'
#' @name CancerRNASig_imports
#' @keywords internal
#' @importFrom stats anova lm predict quantile sd na.omit setNames
#' @importFrom utils data
NULL