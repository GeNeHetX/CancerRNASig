#' PhenoSigDB Integration for CancerRNASig
#'
#' @description 
#' This module provides integration with phenosigdb, allowing CancerRNASig
#' to use the latest curated signatures with clean, normalized naming.
#' 
#' All signatures are returned with human gene names (mouse genes are 
#' automatically translated via phenosigdb's homology mapping).
#'
#' @name phenosigdb_integration
#' @importFrom phenosigdb list_signatures get_signatures
NULL

#' Normalize CancerRNASig signature ID to phenosigdb format
#'
#' @description 
#' Converts legacy CancerRNASig naming conventions to phenosigdb's
#' standardized naming scheme.
#'
#' @param sig_id Character vector of signature IDs in CancerRNASig format
#' @return Character vector of normalized signature IDs in phenosigdb format
#' @keywords internal
.normalize_sig_id <- function(sig_id) {
  if (!is.character(sig_id) || length(sig_id) == 0) {
    return(sig_id)
  }

  # Step 1: Remove .genes suffix
  sig_id <- gsub("\\.genes$", "", sig_id, fixed = FALSE)
  
  # Step 2: Fix known typo
  sig_id <- gsub("Colisson11", "Collisson11", sig_id, fixed = FALSE)
  
  # Step 3: Replace underscores with dots
  sig_id <- gsub("_", ".", sig_id, fixed = FALSE)
  
  # Step 4: Collapse multiple dots
  sig_id <- gsub("\\.+", ".", sig_id, fixed = FALSE)
  
  # Step 5: Prefix normalizations (most specific first)
  sig_id <- gsub("^IMMU\\.MCPcounter\\.Becht16\\.", "IMMUNE.Becht16.", sig_id, fixed = FALSE)
  sig_id <- gsub("^IMMU\\.Tcell\\.Chu23\\.", "IMMUNE.Chu23.", sig_id, fixed = FALSE)
  sig_id <- gsub("^IMMU\\.Neutro\\.Wu24\\.", "IMMUNE.Wu24.", sig_id, fixed = FALSE)
  sig_id <- gsub("^IMMU\\.", "IMMUNE.", sig_id, fixed = FALSE)
  sig_id <- gsub("^IMMUNE\\.monoMac\\.", "IMMUNE.", sig_id, fixed = FALSE)
  sig_id <- gsub("^CAF\\.", "FIBROBLAST.", sig_id, fixed = FALSE)
  sig_id <- gsub("^CCCA\\.", "CANCER_ATLAS.", sig_id, fixed = FALSE)
  sig_id <- gsub("^siNETs\\.", "SINET.", sig_id, fixed = FALSE)
  sig_id <- gsub("^FibroAtlas\\.", "FIBROBLAST.", sig_id, fixed = FALSE)
  sig_id <- gsub("^Organoid\\.Atlas\\.", "ORGANOID.", sig_id, fixed = FALSE)
  sig_id <- gsub("^scGastricTME\\.", "GASTRIC_CANCER.", sig_id, fixed = FALSE)
  sig_id <- gsub("^scIBD\\.", "IBD.", sig_id, fixed = FALSE)
  sig_id <- gsub("^scGC\\.", "GASTRIC_CANCER.", sig_id, fixed = FALSE)
  sig_id <- gsub("^CCK\\.STIM\\.", "CCA.", sig_id, fixed = FALSE)
  sig_id <- gsub("^CCK\\.", "CCA.", sig_id, fixed = FALSE)
  sig_id <- gsub("^Stroma\\.Grunwald21\\.", "PDAC.Grunwald21.", sig_id, fixed = FALSE)
  sig_id <- gsub("^Stroma\\.", "STROMA.", sig_id, fixed = FALSE)
  sig_id <- gsub("^scGI\\.", "GI.", sig_id, fixed = FALSE)
  sig_id <- gsub("^scDuctalPancreas\\.Mouse\\.", "PANCREAS.", sig_id, fixed = FALSE)
  sig_id <- gsub("^scPancreas\\.Mouse\\.", "PANCREAS.", sig_id, fixed = FALSE)
  sig_id <- gsub("^scADM\\.Mouse\\.", "GASTRIC.", sig_id, fixed = FALSE)
  sig_id <- gsub("^SPEM\\.", "GASTRIC.", sig_id, fixed = FALSE)
  sig_id <- gsub("^PDAC\\.PDAssigner\\.", "PDAC.", sig_id, fixed = FALSE)
  sig_id <- gsub("^PDAC\\.PDXph1\\.", "PDAC.", sig_id, fixed = FALSE)

  # Step 6: Handle Greek letters (Unicode to ASCII)
  sig_id <- gsub("\u03b2", "beta", sig_id, fixed = FALSE)  # β
  sig_id <- gsub("\u03b3", "gamma", sig_id, fixed = FALSE)  # γ
  sig_id <- gsub("TGF\u03b2", "TGFbeta", sig_id, fixed = FALSE)
  sig_id <- gsub("IFN\u03b3", "IFNgamma", sig_id, fixed = FALSE)

  # Step 7: Remove trailing dots
  sig_id <- gsub("\\.$", "", sig_id, fixed = FALSE)

  return(sig_id)
}

#' Get signature genes from phenosigdb
#'
#' @description
#' Retrieves gene list for a signature, using CancerRNASig naming convention.
#' Genes are automatically returned in human symbols (mouse genes translated).
#'
#' @param sig_id Signature ID in CancerRNASig format (e.g., "siNETs.Patte25_epithelial.genes")
#' @param reference_species Species for gene names: "human" (default), "mouse", or "original"
#' @return Character vector of gene symbols
#' @export
#' @examples
#' genes <- get_signature_genes("siNETs.Patte25_epithelial.genes")
get_signature_genes <- function(sig_id, reference_species = "human") {
  normalized_id <- .normalize_sig_id(sig_id)
  
  tryCatch({
    sig_list <- phenosigdb::get_signatures(normalized_id, reference_species = reference_species)
    if (is.list(sig_list) && length(sig_list) > 0) {
      return(unlist(sig_list, use.names = FALSE))
    }
    stop("Signature not found in phenosigdb: ", normalized_id)
  }, error = function(e) {
    stop("Failed to get signature '", sig_id, "' (normalized: ", normalized_id, "): ", e$message)
  })
}

#' Get signature metadata from phenosigdb
#'
#' @description
#' Retrieves metadata for a signature, using CancerRNASig naming convention.
#'
#' @param sig_id Signature ID in CancerRNASig format
#' @return Data frame with signature metadata
#' @export
#' @examples
#' meta <- get_signature_metadata("siNETs.Patte25_epithelial.genes")
get_signature_metadata <- function(sig_id) {
  normalized_id <- .normalize_sig_id(sig_id)
  
  tryCatch({
    meta <- phenosigdb::list_signatures(query = normalized_id, reference_species = "human")
    if (!is.null(meta) && nrow(meta) > 0) {
      return(meta)
    }
    stop("Metadata not found in phenosigdb: ", normalized_id)
  }, error = function(e) {
    stop("Failed to get metadata for '", sig_id, "': ", e$message)
  })
}

#' List all CancerRNASig-compatible signatures from phenosigdb
#'
#' @description
#' Returns metadata for all signatures that correspond to CancerRNASig's
#' legacy signatures, loaded from phenosigdb with clean naming.
#'
#' @param reference_species Species for gene names: "human" (default), "mouse", or "original"
#' @return Data frame with signature metadata
#' @export
#' @examples
#' all_sigs <- list_cancer_rna_sig_signatures()
list_cancer_rna_sig_signatures <- function(reference_species = "human") {
  tryCatch({
    # Get all signatures from phenosigdb
    meta <- phenosigdb::list_signatures(reference_species = reference_species)
    
    if (is.null(meta) || nrow(meta) == 0) {
      stop("No signatures found in phenosigdb")
    }
    
    # These are the phenosigdb prefixes that correspond to CancerRNASig signatures
    crna_prefixes <- c(
      "CANCER_ATLAS", "CCA", "FIBROBLAST", "SINET", "ORGANOID",
      "IBD", "GASTRIC_CANCER", "GASTRIC", "IMMUNE", "PDAC",
      "PANCREAS", "STROMA", "GI"
    )
    
    # Filter signatures that match these prefixes
    pattern <- paste0("^(", paste(crna_prefixes, collapse = "|"), ")\\.")
    matching <- grepl(pattern, meta$signature_id)
    
    if (sum(matching) == 0) {
      warning("No signatures matching CancerRNASig prefixes found")
      return(meta[integer(0), ])
    }
    
    return(meta[matching, ])
    
  }, error = function(e) {
    stop("Failed to list signatures: ", e$message)
  })
}
