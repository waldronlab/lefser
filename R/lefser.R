fillPmatZmat <- function(group,
                         block,
                         relab_sub,
                         p.threshold)
{
  if(nrow(relab_sub) == 0L){
    return(relab_sub)
  }
  # creates a list of boolean vectors, each vector indicates
  # existance (TRUE) or absence (FALSE) of a class/sub-class combination
  combos <- apply(
    expand.grid(levels(group), levels(block)), 1L, paste0, collapse = "")
  combined <- paste0(as.character(group), as.character(block))
  logilist <- lapply(setNames(nm = sort(combos)), `==`, combined)

  ## uses Wilcoxon rank-sum test to test for significant differential abundances between
  ## subclasses of one class against subclasses of all othe classes; results are saved in
  ## "pval_mat" and "z_mat" matrices
  whichlist <- lapply(logilist, which)
  sblock <- seq_along(levels(block))
  trelab_sub <- t(relab_sub)
  iters <- expand.grid(sblock, sblock + length(sblock))
  group_formats <- apply(iters, 1L, function(x) {
    ind <- unlist(whichlist[x])
    apply(trelab_sub, 2L, function(g) {
      wx <- suppressWarnings(coin::wilcox_test(g ~ group, subset = ind))
      cbind.data.frame(
          p.value = coin::pvalue(wx), statistic = coin::statistic(wx)
      )
    })
  })

  res <- lapply(group_formats, function(x) do.call(rbind, x))
  pval_mat <- do.call(cbind, lapply(res, `[[`, "p.value"))
  z_mat <- do.call(cbind, lapply(res, `[[`, "statistic"))

  rownames(pval_mat) <- rownames(relab_sub)
  rownames(z_mat) <- rownames(relab_sub)

  ## converts "pval_mat" into boolean matrix "logical_pval_mat" where
  ## p-values <= wilcoxon.threshold
  logical_pval_mat <- pval_mat <= p.threshold * 2.0
  logical_pval_mat[is.na(logical_pval_mat)] <- FALSE

  ## determines which rows (features) have all p-values<=0.05
  ## and selects such rows from the matrix of z-statistics
  sub <- apply(logical_pval_mat, 1L, all)
  z_mat_sub <- z_mat[sub, , drop = FALSE]
    # confirms that z-statistics of a row all have the same sign
  sub <- abs(rowSums(z_mat_sub)) == rowSums(abs(z_mat_sub))
  relab_sub[names(sub[sub]), , drop = FALSE]
}

## ensures that more than half of the values in each for each feature are unique
## if that is not the case then a count value is altered by adding it to a small value
## generated via normal distribution with mean=0 and sd=5% of the count value
createUniqueValues <- function(df, group){
  orderedrows <- rownames(df)
  splitdf <- split(df, group)
  maxim <- vapply(table(group), function(x) max(x * 0.5, 4), numeric(1L))
  for (i in seq_along(splitdf)) {
    sdat <- splitdf[[i]]
    splitdf[[i]][] <- lapply(sdat, function(cols) {
        if (length(unique(cols)) > maxim[i])
            cols
        else
            abs(cols + rnorm(
                length(cols), mean = 0, sd = max(cols * 0.05, 0.01))
            )
    })
  }
  df <- do.call(rbind, unname(splitdf))
  df[match(orderedrows, rownames(df)),, drop = FALSE]
}

contastWithinClassesOrFewPerClass <-
  function(relab_sub_t_df, rand_s, min_cl, ncl, groups) {
    cols <- relab_sub_t_df[rand_s, , drop = FALSE]
    cls <- relab_sub_t_df$class[rand_s]
    # if the number of classes is less than the actual number (typically two)
    # of classes in the dataframe then return TRUE
    if (length(unique(cls)) < ncl) {
      return (TRUE)
    }
    # detect if for each class there are not fewer than the minimum (min_cl) number of samples
    if (any(table(cls) < min_cl)) {
      return (TRUE)
    }
    # separate the randomly selected samples (cols) into a list of the two classes
    drops <- c("class")
    by_class <-
      lapply(seq_along(groups), function(x) {
        cols[cols[, "class"] == groups[x],!(names(cols) %in% drops)]
      })

    # makes sure that within each class all features have at least min_cl unique count values
    for (i in seq_along(groups)) {
      unique_counts_per_microb = apply(by_class[[i]], 2, function(x) {
        length(unique(x))
      })
      if ((any(unique_counts_per_microb <= min_cl) &
           min_cl > 1) |
          (min_cl == 1 & any(unique_counts_per_microb <= 1))) {
        return (TRUE)
      }
    }
    return (FALSE)

  }

ldaFunction <- function (data, lfk, rfk, min_cl, ncl, groups) {
  # test 1000 samples for contrast within classes per feature
  # and that there is at least a minimum number of samples per class
  for (j in 1:1000) {
    rand_s <- sample(seq_len(lfk), rfk, replace = TRUE)
    if (!contastWithinClassesOrFewPerClass(data, rand_s, min_cl, ncl, groups)) {
      break
    }
  }
  # lda with rfk number of samples
  lda.fit <- lda(class ~ ., data = data, subset = rand_s)
  # coefficients that transform observations to discriminants
  w <- lda.fit$scaling[, 1]
  # scaling of lda coefficients
  w.unit <- w / sqrt(sum(w ^ 2))
  sub_d <- data[rand_s,]
  ss <- sub_d[,-match("class", colnames(sub_d))]
  xy.matrix <- as.matrix(ss)
  # the original matrix is transformed
  LD <- xy.matrix %*% w.unit
  # effect size is calculated as difference between averaged disciminants
  # of two classes
  effect_size <-
    abs(mean(LD[sub_d[, "class"] == 1]) - mean(LD[sub_d[, "class"] == 0]))
  # scaling lda coefficients by the efect size
  scal <- w.unit * effect_size
  # mean count values per fclass per feature
  rres <- lda.fit$means

  coeff <- vector("numeric", length(scal))
  for (v in seq_along(scal)) {
    if (!is.na(scal[v])) {
      coeff[v] <- abs(scal[v])
    } else{
      coeff[v] <- 0
    }

  }
  # count value differences between means of two classes for each feature
  lda.means.diff <- (lda.fit$means[2,] - lda.fit$means[1,])
  # difference between a feature's class means and effect size adjusted lda coefficient
  # are averaged for each feature
  (lda.means.diff + coeff) / 2
}

filterKruskal <- function(relab, group, p.value) {
  # applies "kruskal.test.alt" function to each row (feature) of relab
  # to detect differential abundance between classes, 0 and 1
  kw.res <- apply(relab, 1L, function(x) {
    kruskal.test(x ~ group)[["p.value"]]
  })
  # TRUE for p-values less than or equal to kw.threshold
  kw.sub <- kw.res < p.value

  # NAs are FALSE
  kw.sub[is.na(kw.sub)] <- FALSE

  # extracts features with statistically significant differential abundance
  # from "relab" matrix
  relab[kw.sub,]
}

#' R implementation of the LEfSe method
#'
#' Perform a LEfSe analysis: the function carries out differential analysis
#' between two sample groups for multiple microorganisms and uses linear discirminant analysis
#' to establish their effect sizes. Subclass information for each class can be incorporated
#' into the analysis (see examples). Microorganisms with large differences between two sample groups
#' are identified as biomarkers.
#'
#' @details
#' The LEfSe method expects relative abundances in the `expr` input. A warning
#' will be emitted if the column sums do not result in 1. Use the `relativeAb`
#' helper function to convert the data in the `SummarizedExperiment` to relative
#' abundances. The `checkAbundances` argument enables checking the data
#' for presence of relative abundances and can be turned off by setting the
#' argument to `FALSE`.
#'
#' @param relab A [SummarizedExperiment-class] with relative
#'   abundances in the assay
#' @param expr (`deprecated`) Use `relab` instead. A [SummarizedExperiment-class]
#'   with relative abundances in the assay
#' @param kruskal.threshold numeric(1) The p-value for the Kruskal-Wallis Rank
#' Sum Test (default 0.05).
#' @param wilcox.threshold numeric(1) The p-value for the Wilcoxon Rank-Sum Test
#' when 'blockCol' is present (default 0.05).
#' @param lda.threshold numeric(1) The effect size threshold (default 2.0).
#' @param groupCol character(1) Column name in `colData(relab)` indicating
#' groups, usually a factor with two levels (e.g., `c("cases", "controls")`;
#' default "GROUP").
#' @param blockCol character(1) Optional column name in `colData(relab)`
#' indicating the blocks, usually a factor with two levels (e.g.,
#' `c("adult", "senior")`; default NULL).
#' @param assay The i-th assay matrix in the `SummarizedExperiment` ('relab';
#' default 1).
#' @param trim.names If `TRUE` extracts the most specific taxonomic rank of organism.
#' @param checkAbundances `logical(1)` Whether to check if the assay data in the
#'   `relab` input are relative abundances or counts. If counts are found, a
#'   warning will be emitted (default `TRUE`).
#' @param \ldots Additional inputs to lower level functions (not used).
#' @return
#' The function returns a `data.frame` with two columns, which are
#' names of microorganisms and their LDA scores.
#'
#' @importFrom stats kruskal.test reorder rnorm
#' @importFrom coin pvalue statistic wilcox_test
#' @importFrom MASS lda
#' @importFrom methods as is
#' @importFrom stats setNames
#' @importFrom utils tail
#' @import SummarizedExperiment
#'
#' @examples
#'     # (1) Using classes only
#'     data(zeller14)
#'     # exclude 'adenoma'
#'     zeller14 <- zeller14[, zeller14$study_condition != "adenoma"]
#'     res_group <- lefser(zeller14, groupCol = "study_condition")
#'     head(res_group)
#'
#'     # (2) Using classes and sublasses
#'     data(zeller14)
#'     # exclude 'adenoma'
#'     zeller14 <- zeller14[, zeller14$study_condition != "adenoma"]
#'     res_block <- lefser(
#'          zeller14, groupCol = "study_condition", blockCol = "age_category"
#'     )
#'     head(res_block)
#' @export
lefser <-
  function(relab,
           kruskal.threshold = 0.05,
           wilcox.threshold = 0.05,
           lda.threshold = 2.0,
           groupCol = "GROUP",
           blockCol = NULL,
           assay = 1L,
           trim.names = FALSE,
           checkAbundances = TRUE,
           ...,
           expr
) {
    if (!missing(expr)) {
        .Deprecated(
            msg = "The 'expr' argument is deprecated, use 'relab' instead."
        )
        relab_data <- assay(expr, i = assay)
    } else {
        relab_data <- assay(relab, i = assay)
    }
    if (checkAbundances && !identical(all.equal(colSums(relab_data), rep(1e6, ncol(relab_data)), check.attributes = FALSE), TRUE)) {
        stop("Convert counts to relative abundances with 'relativeAb()'")
    }
        
    groupf <- colData(relab)[[groupCol]]
    groupf <- as.factor(groupf)
    lgroupf <- levels(groupf)
    if (is.null(groupf) || !identical(length(lgroupf), 2L)) {
        stop(
            "'groupCol' must refer to a valid dichotomous (two-level) variable"
        )
    }
    message(
        "The outcome variable is specified as '", groupCol,
        "' and the reference category is '", lgroupf[1],
        "'.\n See `?factor` or `?relevel` to change the reference category."
    )
    relab_sub <- filterKruskal(relab = relab_data, group = groupf, p.value = kruskal.threshold)

    if (!is.null(blockCol)) {
        block <- as.factor(colData(relab)[[blockCol]])
        block <- droplevels(block)
        relab_sub <- fillPmatZmat(group = groupf, block = block, relab_sub = relab_sub, p.threshold = wilcox.threshold)
    }
    
    if(nrow(relab_sub) == 0L){
      return(.return_no_results())
    }

    # transposes matrix and add a "class" (i.e., groupf) column
    # matrix converted to dataframe
    relab_sub_t <- t(relab_sub)
    relab_sub_t_df <- as.data.frame(relab_sub_t)
    relab_sub_t_df <- createUniqueValues(df = relab_sub_t_df, group = groupf)
    relab_sub_t_df <- cbind(relab_sub_t_df, class = groupf)

    # number of samples (i.e., subjects) in the dataframe
    lfk <- nrow(relab_sub_t_df)
    # rfk is the number of subject that will be used in linear discriminant analysis
    rfk <- floor(lfk * 2 / 3)
    # number of classes (two-levels)
    ncl <- length(lgroupf)
    # count samples in each class of the dataframe, select the number from the class with a smaller
    # count of samples and multiply that number by 2/*2/3*0.5
    min_cl <-
      as.integer(min(table(relab_sub_t_df$class)) * 2 / 3 * 2 / 3 *
                   0.5)
    # if min_cl is less than 1, then make it equal to 1
    min_cl <- max(min_cl, 1)

    # lda_fn repeated 30 times, producing a matrix of 30 scores per feature
    eff_size_mat <-
      replicate(30, suppressWarnings(ldaFunction(
        relab_sub_t_df, lfk, rfk, min_cl, ncl, lgroupf
      )), simplify = TRUE)

    # mean of 30 scores per feature
    raw_lda_scores <- rowMeans(eff_size_mat)

    # processing of score
    processed_scores <-
      sign(raw_lda_scores) * log((1 + abs(raw_lda_scores)), 10)

    # sorting of scores
    processed_sorted_scores <- sort(processed_scores)
    scores_df <- data.frame(Names = names(processed_sorted_scores),
                            scores = as.vector(processed_sorted_scores),
                            stringsAsFactors = FALSE)

    scores_df <- .trunc(scores_df, trim.names)

    threshold_scores <- abs(scores_df$scores) >= lda.threshold
    res_scores <- scores_df[threshold_scores, ]
    class(res_scores) <- c("lefser_df", class(res_scores))
    attr(res_scores, "groups") <- lgroupf
    if(nrow(res_scores) == 0L){
      return(.return_no_results())
    }
    res_scores
  }

.return_no_results <- function() {
  message("No significant features found.")
  res_scores <- data.frame(Names=character(), scores=numeric())
  class(res_scores) <- c("lefser_df", class(res_scores))
  return(res_scores)
}