#' Title
#' Perform a LEfSe analysis
#'
#' @param expr
#' Explain expr
#'
#' @return
#' Explain what the function returns
#'
#' @export
#' @importFrom stats kruskal.test reorder rnorm
#' @importFrom coin pvalue statistic wilcox_test
#' @importFrom MASS lda
#' @importFrom methods as is
#' @import SummarizedExperiment

#'
#' @param expr
#' @param kw.threshold
#' @param wilcoxon.threshold
#'
#' @examples
#' data(zeller14)
#' zeller14 <- zeller14[, zeller14$study_condition != "adenoma"]
#' zeller14$GROUP <- ifelse(zeller14$study_condition == "control", 0, 1)
#' # keep only taxa with non-zero counts in at least 10 specimens:
#' results <- lefserAnalysis(zeller14)
#' head(results)

lefserAnalysis <- function (expr, kw.threshold = 0.05, wilcoxon.threshold = 0.05)

{
  if (is(expr, "ExpressionSet"))
    expr <- as(expr, "SummarizedExperiment")
  isSE <- is(expr, "SummarizedExperiment")
  if (isSE) {
    se <- expr
    info <- EnrichmentBrowser:::.extractInfoFromSE(se)
    expr <- info$expr
    grp <- info$grp
    blk <- info$blk
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
  group <- factor(grp)
  block <- factor(blk)

  # extracts p-values from Kruskal-Wallis Rank Sum Test
  kruskal.test.alt <- function(x) {
    kruskal.test(x ~ group)$p.value
  }

  # applies "kruskal.test.alt" function to each row (feature) of expr
  # to detect differential abundance between classes, 0 and 1
  kw.res <- apply(expr, 1, kruskal.test.alt)

  # selects p-values less than or equal to kw.threshold
  kw.sub <- kw.res <= kw.threshold

  # eliminates NAs
  kw.sub[is.na(kw.sub)] <- FALSE

  # extracts features with statistically significant differential abundance
  # from "expr" matrix
  expr1 <- expr[kw.sub, ]

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
    matrix(0, nrow(expr1), length(unique(block)) ^ length(unique(group)))
  z_mat <-
    matrix(0, nrow(expr1), length(unique(block)) ^ length(unique(group)))

  rownames(pval_mat) <- rownames(expr1)
  rownames(z_mat) <- rownames(expr1)

  # uses Wilcoxon rank-sum test to test for significant differential abundances between
  # subclasses of one class against subclasses of all othe classes; results are saved in
  # "pval_mat" and "z_mat" matrices
  c <- 0
  for (i in seq_along(unique(block))) {
    for (j in seq_along(unique(block))) {
      j <- j + length(unique(block))
      mat = cbind(expr1[, logical.list.of.subclasses[[i]]],
                  expr1[, logical.list.of.subclasses[[j]]])
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
  expr1_sub <- expr1[names(sub[sub == TRUE]),]

  # transposes matrix and add a "class" (i.e., group) column
  # matrix converted to dataframe
  expr1_sub_t <- t(expr1_sub)
  expr1_sub_t <- cbind(expr1_sub_t, class = (as.numeric(group) - 1))
  expr1_sub_t_df <- data.frame(expr1_sub_t)

  # makes sure that more than half of the values in each for each feature are unique
  # if that is not the case then a count value is altered by adding it to a small value
  # generated via normal distribution with mean=0 and sd=5% of the count value
  df <- expr1_sub_t_df
  for (i in seq_along(c(1:(ncol(df) - 1)))) {
    for (j in seq_along(groups)) {
      equality = df[df[, "class"] == groups[j], i]
      if (length(unique(equality)) > max(length(equality) * 0.5, 4)) {
        next
      } else{
        for (k in seq_along(equality)) {
          equality[k] = abs(equality[k] + rnorm(
            1,
            mean = 0,
            sd = max(equality[k] * 0.05, 0.01)
          ))
        }

      }
      df[df[, "class"] == groups[j], i] = equality
    }
  }

  expr1_sub_t_df <- df
  # number of samples (i.e., subjects) in the dataframe
  lfk <- nrow(expr1_sub_t_df)
  # rfk is the number of subject that will be used in linear discriminant analysis
  rfk <- floor(lfk * 2 / 3)
  # number of classes (typically two)
  ncl <- length(groups)
  # count samples in each class of the dataframe, select the number from the class with a smaller
  # count of samples and multiply that number by 2/*2/3*0.5
  min_cl <- as.integer(min(table(expr1_sub_t_df$class)) * 2 / 3 * 2 / 3 *
                         0.5)
  # if min_cl is less than 1, then make it equal to 1
  min_cl <- max(min_cl, 1)


  contast_within_classes_or_few_per_class <-
    function(expr1_sub_t_df, rand_s, min_cl, ncl) {
      cols <- expr1_sub_t_df[rand_s, ]
      cls <- expr1_sub_t_df$class[rand_s]
      # if the number of classes is less than the actual number (typically two)
      # of classes in the dataframe then return TRUE
      if (length(unique(cls)) < ncl) {
        return (TRUE)
      }
      # detect if for each class there are not fewer than the minimum (min_cl) number of samples
      if (TRUE %in% c(table(cls) < min_cl)) {
        return (TRUE)
      }
      # separate the randomly selected samples (cols) into a list of the two classes
      drops <- c("class")
      by_class <-
        lapply(seq_along(groups), function(x) {
          cols[cols[, "class"] == groups[x], !(names(cols) %in% drops)]
        })

      # makes sure that within each class all features have at least min_cl unique count values
      for (i in seq_along(groups)) {
        unique_counts_per_microb = apply(by_class[[i]], 2, function(x) {
          length(unique(x))
        })
        if ((TRUE %in% c(unique_counts_per_microb <= min_cl) &
             min_cl > 1) |
            (min_cl == 1 & (TRUE %in% c(unique_counts_per_microb <= 1)))) {
          return (TRUE)
        }
      }
      return (FALSE)

    }

  lda_fn <- function (data) {
    # test 1000 samples for contrast within classes per feature
    # and that there is at least a minimum number of samples per class
    for (j in seq_along(1:1000)) {
      rand_s <- sample(c(1:lfk), rfk, replace = TRUE)
      if (!contast_within_classes_or_few_per_class(data, rand_s, min_cl, ncl)) {
        break
      }
    }
    # lda with rfk number of samples
    lda.fit <- lda(class ~ ., data = data, subset = rand_s)
    # coefficients that transform observations to discriminants
    w <- lda.fit$scaling[, 1]
    # scaling of lda coefficients
    w.unit <- w / sqrt(sum(w ^ 2))
    sub_d <- data[rand_s, ]
    ss <- sub_d[, -match("class", colnames(sub_d))]
    xy.matrix <- as.matrix(ss)
    # the original matrix is transformed
    LD <- xy.matrix %*% w.unit
    # effect size is calculated as difference between averaged disciminants
    # of two classes
    effect_size <-
      abs(mean(LD[sub_d[, "class"] == 0]) - mean(LD[sub_d[, "class"] == 1]))
    # scaling lda coefficients by the efect size
    scal <- w.unit * effect_size
    # mean count values per fclass per feature
    rres <- lda.fit$means
    rowns <- rownames(rres)
    lenc <- length(colnames(rres))

    coeff <- vector("numeric", length(scal))
    for (v in seq_along(scal)) {
      if (is.na(scal[v]) != TRUE) {
        coeff[v] <- abs(scal[v])
      } else{
        coeff[v] <- 0
      }

    }
    # count value differences between means of two classes for each feature
    lda.means.diff <- (lda.fit$means[2, ] - lda.fit$means[1, ])
    # difference between a feature's class means and effect size adjusted lda coefficient
    # are averaged for each feature
    (lda.means.diff + coeff) / 2
  }

  # lda_fn repeated 30 times, producing a matrix of 30 scores per feature
  eff_size_mat <- replicate(30, lda_fn(expr1_sub_t_df), simplify = T)

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
  return(scores_df)
}
