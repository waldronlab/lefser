#' Utility function to calculate relative abundances
#'
#' @description
#' The function calculates the column totals and divides each value within the
#' column by the respective column total.
#'
#' @inheritParams lefser
#'
#' @param se A SummarizedExperiment object with counts
#'
#' @examples
#'
#' se <- SummarizedExperiment(
#'     assays = list(
#'         counts = matrix(
#'             rep(1, 4), ncol = 1, dimnames = list(LETTERS[1:4], "SAMP")
#'         )
#'     )
#' )
#' assay(se)
#' assay(relativeAb(se))
#'
#' @description
#' This function calculates the relative abundance of each feature in the SummarizedExperiment 
#' object containing count data, expressed as counts per million (CPM) 
#' 
#' @returns returns a new SummarizedExperiment object with counts per million
#' calculated and added as a new assay named rel_abs.
#' 
#' @export
relativeAb <- function(se, assay = 1L) {
  assay_data <- assay(se, i = assay)
  csums <- colSums(assay_data)
  div <- matrix(rep(csums, each = nrow(assay_data)), ncol = ncol(assay_data))
  res <- assay_data / div * 1e6
  assaylist <- assays(se)
  newalist <- append(
    assaylist, values = S4Vectors::SimpleList(rel_abs = res), after = 0L
  )
  assays(se) <- newalist
  se
}

# Truncate the feature name
.trunc <- function(scores_df, trim.names){
    Names <- gsub("`", "", scores_df[["Names"]])
    if (trim.names) {
        listNames <- strsplit(Names, "\\||\\.")
        Names <- vapply(listNames, tail, character(1L), 1L)
    }
    scores_df[["Names"]] <- Names
    return(scores_df)
}
