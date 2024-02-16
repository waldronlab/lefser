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
#' @param assay Name or number of the assay of se containing count data to be used for the relative abundance calculation (default: 1L).
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
