#' @inheritParams lefser
#' @param se A SummarizedExperiment object with counts
#' @export
relativeAb <- function(se, assay = 1L) {
    expr_data <- assay(se, i = assay)
    csums <- colSums(expr_data)
    div <- matrix(rep(csums, each = nrow(expr_data)), ncol = ncol(expr_data))
    res <- expr_data / div
    assay(se, i = assay) <- res
    se
}
