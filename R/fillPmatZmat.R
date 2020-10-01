fillPmatZmat <- function(group, block, expr_sub, wilcoxon.threshold)
{
  # extracts p-values from Wilcoxon rank-sum test
  wilcox_test_pvalue <- function(x, group) {
    pvalue(wilcox_test(x ~ group))
  }
  
  # extracts z-statistic from Wilcoxon rank-sum test
  wilcox_test_z <- function(x, group) {
    statistic(wilcox_test(x ~ group))
  }
  # creates a list of boolean vectors, each vector indicates
  # existance (TRUE) or absence (FALSE) of a class/sub-class combination
  logical.list.of.subclasses <- list()
  for (i in seq_along(unique(group))) {
    for (j in seq_along(unique(block))) {
      logical.list.of.subclasses <- append(logical.list.of.subclasses,
                                           list(group == sort(unique(group))[i] &
                                                  block == sort(unique(block))[j]))
    }
  }
  
  # creates empty matrices that will be filled by p-values and z-statistics from Wilcoxon rank-sum test
  pval_mat <-
    matrix(0, nrow(expr_sub), length(unique(block)) ^ length(unique(group)))
  z_mat <-
    matrix(0, nrow(expr_sub), length(unique(block)) ^ length(unique(group)))
  
  rownames(pval_mat) <- rownames(expr_sub)
  rownames(z_mat) <- rownames(expr_sub)

  # uses Wilcoxon rank-sum test to test for significant differential abundances between
  # subclasses of one class against subclasses of all othe classes; results are saved in
  # "pval_mat" and "z_mat" matrices
  c <- 0
  for (i in seq_along(unique(block))) {
    for (j in seq_along(unique(block))) {
      j <- j + length(unique(block))
      mat = cbind(expr_sub[, logical.list.of.subclasses[[i]]],
                  expr_sub[, logical.list.of.subclasses[[j]]])
      group_for_mat = factor(c(rep(
        0, count(logical.list.of.subclasses[[i]] == TRUE)
      ),
      rep(
        1, count(logical.list.of.subclasses[[j]] == TRUE)
      )))
      pval_mat[, c + 1] <-
        apply(mat, 1, wilcox_test_pvalue, group = group_for_mat)
      z_mat[, c + 1] <-
        apply(mat, 1, wilcox_test_z, group = group_for_mat)
      c = c + 1
    }
  }
  
  # number of p-values per feature
  num_comp <- length(unique(block)) ^ length(unique(group))
  
  # converts "pval_mat" into boolean matrix "logical_pval_mat" where p-values <= wilcoxon.threshold
  logical_pval_mat <- pval_mat <= wilcoxon.threshold
  
  # determines which rows (features) have all p-values<=0.05
  # and selects such rows from the matrix of z-statistics
  sub <- which(rowSums(logical_pval_mat) == num_comp)
  z_mat_sub <- z_mat[sub, ]
  
  # confirms that z-statistics of a row all have the same sign
  sub <- abs(rowSums(z_mat_sub)) == abs(rowSums(abs(z_mat_sub)))
  expr_sub <- expr_sub[names(sub[sub == TRUE]),]
}
