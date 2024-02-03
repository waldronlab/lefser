#' Utility function to calculate relative abundances
#'
#' @description
#' The function calculates the column totals and divides each value within the
#' column by the respective column total.
#'
#' @param se A SummarizedExperiment object with counts
#'
#' @return A modified SummarizedExperiment object
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
#' @export
relativeAb <- function(se, assay = 1L) {
  
  assay_data <- assay(se, i = assay)
  
  
  csums <- colSums(assay_data)
  
  
  cat("")
  print(assay_data)
  
  
  res <- assay_data / csums * 1e6
  
  
  cat("")
  print(res)
  
  
  assaylist <- assays(se)
  for (i in seq_along(assaylist)) {
    assaylist[[i]] <- assaylist[[i]] * 1e6
  }
  
  
  newalist <- append(
    assaylist, 
    values = S4Vectors::SimpleList(rel_abs = res), 
    after = 0L
  )
  
  
  assays(se) <- newalist
  
  
  cat("")
  print(se)
  
  
  return(se)
}


se <- SummarizedExperiment(
  assays = list(
    counts = matrix(
      rep(1, 4), ncol = 1, dimnames = list(LETTERS[1:4], "SAMP")
    )
  )
)

assay(se)
assay(relativeAb(se))


