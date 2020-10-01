#' Title
#' Perform a LEfSe analysis
#'
#' @param expr
#' The expr is an \code{\linkS4class{ExpressionSet}} or a \code{\linkS4class{SummarizedExperiment}}.
#' GROUP column should be assigned to meta-data of
#' ExpressionSet/SummarizedExperiment
#' as a class in pData/colData before using lefserAnalysis function.
#' Use '0' and '1' for unaffected (controls) and
#' affected (cases) samples, respectively. Optionally, any number of subclasses can
#' be defined by adding BLOCK column to meta-data of ExpressionSet/SummarizedExperiment.
#' @param kw.threshold
#' The p-value threshold for Kruskal-Wallis Rank Sum Test.
#' The default is at <= 0.05.
#' @param wilcoxon.threshold
#' The p-value for Wilcoxon Rank-Sum Test.
#' The default is at <= 0.05.
#' @param lda.threshold
#' The effect size threshold.
#' The default is at 2.0.
#' @return
#' The function returns a dataframe with two columns, which are
#' names of microorganisms and their LDA scores.
#'
#' @export
#' @importFrom stats kruskal.test reorder rnorm
#' @importFrom coin pvalue statistic wilcox_test
#' @importFrom MASS lda
#' @importFrom methods as is
#' @import SummarizedExperiment
#'
#' @examples
#'     # (1) Using classes only
#'     data(zeller14)
#'     # exclude 'adenoma'
#'     zeller14 <- zeller14[, zeller14$study_condition != "adenoma"]
#'     # assign '0' class to 'conrol' and '1' to 'CRC' (i.e., colorectal cancer)
#'     zeller14$GROUP <- ifelse(zeller14$study_condition == "control", 0, 1)
#'     results <- lefserAnalysis(zeller14)
#'     head(results)
#'
#'     # (2) Using classes and sublasses
#'     data(zeller14)
#'     # exclude 'adenoma'
#'     zeller14 <- zeller14[, zeller14$study_condition != "adenoma"]
#'     # assign '0' class to 'conrol' and '1' to 'CRC' (i.e., colorectal cancer)
#'     zeller14$GROUP <- ifelse(zeller14$study_condition == "control", 0, 1)
#'     # assign '0' class to 'adult' and '1' to 'senior'
#'     zeller14$BLOCK <- ifelse(zeller14$age_category == "adult", 0, 1)
#'     results <- lefserAnalysis(zeller14)
#'     head(results)

lefserAnalysis <- function (expr, kw.threshold = 0.05, wilcoxon.threshold = 0.05, lda.threshold = 2.0)

{
  groupBlockMatrix <- createGroupBlock(expr)
  group <- groupBlockMatrix$group
  block <- groupBlockMatrix$block
  expr <- groupBlockMatrix$expr

  # extracts p-values from Kruskal-Wallis Rank Sum Test
  kruskal.test.alt <- function(x, group) {
    kruskal.test(x ~ group)$p.value
  }

  # applies "kruskal.test.alt" function to each row (feature) of expr
  # to detect differential abundance between classes, 0 and 1
  kw.res <- apply(expr, 1, kruskal.test.alt, group=group)

  # selects p-values less than or equal to kw.threshold
  kw.sub <- kw.res <= kw.threshold

  # eliminates NAs
  kw.sub[is.na(kw.sub)] <- FALSE

  # extracts features with statistically significant differential abundance
  # from "expr" matrix
  expr_sub <- expr[kw.sub, ]

 if (length(block)!=0){
   expr_sub <- fillPmatZmat(group, block, expr_sub, wilcoxon.threshold)
 }

  # transposes matrix and add a "class" (i.e., group) column
  # matrix converted to dataframe
  expr_sub_t <- t(expr_sub)
  expr_sub_t <- cbind(expr_sub_t, class = (as.numeric(group) - 1))
  expr_sub_t_df <- data.frame(expr_sub_t)

  expr_sub_t_df <- createUniqueValues(expr_sub_t_df, groups)

  # number of samples (i.e., subjects) in the dataframe
  lfk <- nrow(expr_sub_t_df)
  # rfk is the number of subject that will be used in linear discriminant analysis
  rfk <- floor(lfk * 2 / 3)
  # number of classes (typically two)
  ncl <- length(groups)
  # count samples in each class of the dataframe, select the number from the class with a smaller
  # count of samples and multiply that number by 2/*2/3*0.5
  min_cl <- as.integer(min(table(expr_sub_t_df$class)) * 2 / 3 * 2 / 3 *
                         0.5)
  # if min_cl is less than 1, then make it equal to 1
  min_cl <- max(min_cl, 1)

  # lda_fn repeated 30 times, producing a matrix of 30 scores per feature
  eff_size_mat <- replicate(30, suppressWarnings(ldaFunction(expr_sub_t_df, lfk, rfk, min_cl, ncl)), simplify = T)

  # mean of 30 scores per feature
  raw_lda_scores <- rowMeans(eff_size_mat)

  # processing of score
  processed_scores <-
    sign(raw_lda_scores) * log((1 + abs(raw_lda_scores)), 10)

  # sorting of scores
  processed_sorted_scores <- sort(processed_scores)

  # extracting the most specific taxonomic rank of an organism
  Names <- vector("character", length(processed_sorted_scores))
  for (i in seq_along(processed_sorted_scores)) {
    vec_of_strings = unlist(strsplit(names(processed_sorted_scores)[i], "[.]"))
    Names[i] = vec_of_strings[length(vec_of_strings)]
  }

  # return a dataframe of taxon name and effect size score
  scores <- as.vector(processed_sorted_scores)
  scores_df <- data.frame(Names, scores)
  threshold_scores <- abs(scores_df$scores)>=lda.threshold
  scores_df <- scores_df[threshold_scores,]
  rownames(scores_df) <- c(1:nrow(scores_df))
  return(scores_df)
}
