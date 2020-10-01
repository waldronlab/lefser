createGroupBlockGroups <- function(expr){
  
  if (is(expr, "ExpressionSet")){
    expr <- as(expr, "SummarizedExperiment")
  }
  isSE <- is(expr, "SummarizedExperiment")
  if (isSE) {
    se <- expr
    expr <- assay(se)
    grp <- colData(se)$GROUP
    blk <- colData(se)$BLOCK
  }
  if (!is.matrix(expr))
    stop(
      paste(
        "Expression data in 'expr' must be either",
        "a matrix, a SummarizedExperiment, or an ExpressionSet"
      )
    )
  if (is.null(grp))
    stop("Group assignment 'grp' must be specified")
  groups <- sort(unique(grp))
  if (!all(groups == c(0, 1)))
    stop(
      paste0(
        "Group classification is not binary:\n",
        "Expected (0, 1) but found (",
        paste(groups, collapse = ", "),
        ")"
      )
    )
  list(group = factor(grp),
         block = factor(blk), expr = expr, groups=groups)
}
