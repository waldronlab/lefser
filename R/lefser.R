fillPmatZmat <- function(fun,
                         group,
                         block,
                         expr_sub,
                         p.threshold)
{
  hasBlocks <- length(levels(block)) > 1L
  # creates a list of boolean vectors, each vector indicates
  # existance (TRUE) or absence (FALSE) of a class/sub-class combination
  combos <- apply(
    expand.grid(levels(group), levels(block)), 1L, paste0, collapse = "")
  combined <- paste0(as.character(group), as.character(block))
  logilist <- lapply(setNames(nm = sort(combos)), `==`, combined)

  pv <- if (hasBlocks) { coin::pvalue } else function(x) { x$p.value }
  stat <- if (hasBlocks) { coin::statistic } else function(x) { x$statistic }
  ## uses Wilcoxon rank-sum test to test for significant differential abundances between
  ## subclasses of one class against subclasses of all othe classes; results are saved in
  ## "pval_mat" and "z_mat" matrices
  whichlist <- lapply(logilist, which)
  sblock <- seq_along(levels(block))
  texp_sub <- t(expr_sub)
  iters <- expand.grid(sblock, sblock + length(sblock))
  group_formats <- apply(iters, 1L, function(x) {
    ind <- unlist(whichlist[x])
    apply(texp_sub, 2L, function(g) {
      wx <- suppressWarnings(fun(g ~ group, subset = ind))
      cbind.data.frame(p.value = pv(wx), statistic = stat(wx))
    })
  })

  res <- lapply(group_formats, function(x) do.call(rbind, x))
  pval_mat <- do.call(cbind, lapply(res, `[[`, "p.value"))
  z_mat <- do.call(cbind, lapply(res, `[[`, "statistic"))

  rownames(pval_mat) <- rownames(expr_sub)
  rownames(z_mat) <- rownames(expr_sub)

  ## converts "pval_mat" into boolean matrix "logical_pval_mat" where
  ## p-values <= wilcoxon.threshold
  logical_pval_mat <- pval_mat <= p.threshold
  logical_pval_mat[is.na(logical_pval_mat)] <- FALSE

  ## determines which rows (features) have all p-values<=0.05
  ## and selects such rows from the matrix of z-statistics
  sub <- apply(logical_pval_mat, 1L, all)
  if (hasBlocks) {
    z_mat_sub <- z_mat[sub,]
    # confirms that z-statistics of a row all have the same sign
    sub <- abs(rowSums(z_mat_sub)) == rowSums(abs(z_mat_sub))
  }
  expr_sub[names(sub[sub]), ]
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
  df[match(orderedrows, rownames(df)),]
}

contastWithinClassesOrFewPerClass <-
  function(expr_sub_t_df, rand_s, min_cl, ncl, groups) {
    cols <- expr_sub_t_df[rand_s,]
    cls <- expr_sub_t_df$class[rand_s]
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
    abs(mean(LD[sub_d[, "class"] == 0]) - mean(LD[sub_d[, "class"] == 1]))
  # scaling lda coefficients by the efect size
  scal <- w.unit * effect_size
  # mean count values per fclass per feature
  rres <- lda.fit$means
  rowns <- rownames(rres)
  lenc <- length(colnames(rres))

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

.numeric01 <- function(x) {
    x <- as.factor(x)
    uvals <- levels(x)
    ifelse(x == uvals[1L], 0L, 1L)
}

#' R implementation of the LEfSe method
#'
#' Perform a LEfSe analysis: the function carries out differential analysis
#' between two sample groups for multiple microorganisms and uses linear discirminant analysis
#' to establish their effect sizes. Subclass information for each class can be incorporated
#' into the analysis (see examples). Microorganisms with large differences between two sample groups
#' are identified as biomarkers.
#'
#' @param expr A \code{\linkS4class{SummarizedExperiment}} with expression data.
#' @param groupCol character(1) Column name in `colData(expr)` indicating
#' groups, usually a factor with two levels (e.g., `c("cases", "controls")`;
#' default "GROUP").
#' @param blockCol character(1) Column name in `colData(expr)` indicating the
#' blocks, usually a factor with two levels (e.g., `c("adult", "senior")`;
#' default NULL).
#' @param assay The i-th assay matrix in the `SummarizedExperiment` ('expr').
#' Defaults to 1.
#' @param p.threshold numeric(1) The p-value for the Kruskal-Wallis Rank Sum
#' Test. If 'blockCol' is provided, it indicates the p-value for the Wilcoxon
#' Rank-Sum Test (default 0.05).
#' @param lda.threshold
#' The effect size threshold.
#' The default is at 2.0.
#' @return
#' The function returns a dataframe with two columns, which are
#' names of microorganisms and their LDA scores.
#'
#' @importFrom stats kruskal.test reorder rnorm
#' @importFrom coin pvalue statistic wilcox_test
#' @importFrom MASS lda
#' @importFrom methods as is
#' @importFrom stats setNames
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
  function(expr,
           p.threshold = 0.05,
           lda.threshold = 2.0,
           groupCol = "GROUP",
           blockCol = NULL,
           assay = 1L)
  {
    groupf <- colData(expr)[[groupCol]]
    if (is.null(groupf))
        stop("A valid group assignment 'groupCol' must be provided")
    groupf <- as.factor(groupf)
    groupsf <- levels(groupf)
    if (length(groupsf) != 2L)
      stop(
        "Group classification is not dichotomous:\n",
        "Found (", paste(groupsf, collapse = ", "), ")"
      )
    group <- .numeric01(groupf)
    groups <- 0:1

    block <- factor(rep(1L, length(group)))
    fun <- stats::kruskal.test
    if (!is.null(blockCol)) {
        block <- as.factor(colData(expr)[[blockCol]])
        fun <- coin::wilcox_test
    }

    expr <- assay(expr, i = assay)

    expr_sub <- fillPmatZmat(fun, groupf, block, expr, p.threshold)

    # transposes matrix and add a "class" (i.e., group) column
    # matrix converted to dataframe
    expr_sub_t <- t(expr_sub)
    expr_sub_t_df <- as.data.frame(expr_sub_t)
    expr_sub_t_df <- createUniqueValues(expr_sub_t_df, groupf)
    expr_sub_t_df <- cbind(expr_sub_t_df, class = group)

    # number of samples (i.e., subjects) in the dataframe
    lfk <- nrow(expr_sub_t_df)
    # rfk is the number of subject that will be used in linear discriminant analysis
    rfk <- floor(lfk * 2 / 3)
    # number of classes (two)
    ncl <- length(groups)
    # count samples in each class of the dataframe, select the number from the class with a smaller
    # count of samples and multiply that number by 2/*2/3*0.5
    min_cl <-
      as.integer(min(table(expr_sub_t_df$class)) * 2 / 3 * 2 / 3 *
                   0.5)
    # if min_cl is less than 1, then make it equal to 1
    min_cl <- max(min_cl, 1)

    # lda_fn repeated 30 times, producing a matrix of 30 scores per feature
    eff_size_mat <-
      replicate(30, suppressWarnings(ldaFunction(
        expr_sub_t_df, lfk, rfk, min_cl, ncl, groups
      )), simplify = TRUE)

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
    threshold_scores <- abs(scores_df$scores) >= lda.threshold
    scores_df <- scores_df[threshold_scores, ]
    rownames(scores_df) <- seq_len(nrow(scores_df))
    return(scores_df)
  }
