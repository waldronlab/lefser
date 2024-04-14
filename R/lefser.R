## Wilcoxon Rank-Sum Test for the sub-classes
fillPmatZmat <- function(group,
                         block,
                         relab_sub,
                         p.threshold)
{
  if(nrow(relab_sub) == 0L){
    return(relab_sub)
  }
  # creates a list of boolean vectors, each vector indicates
  # existence (TRUE) or absence (FALSE) of a class/sub-class combination
  combos <- apply(
    expand.grid(levels(group), levels(block)), 1L, paste0, collapse = "")
  combined <- paste0(as.character(group), as.character(block))
  logilist <- lapply(setNames(nm = sort(combos)), `==`, combined)

  ## uses Wilcoxon rank-sum test to test for significant differential abundances between
  ## subclasses of one class against subclasses of all other classes; results are saved in
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

## Quality control of the random subset for LDA
## This function returns `FALSE` when the selected random samples meet 
## the 3 quality criteria
contastWithinClassesOrFewPerClass <- function(relab_sub_t_df, 
                                              rand_s, 
                                              min_cl, 
                                              ncl, 
                                              groups) {
    cols <- relab_sub_t_df[rand_s, , drop = FALSE]
    cls <- relab_sub_t_df$class[rand_s]
    
    ## 1. Minimum number of classes in the random subset
    if (length(unique(cls)) != ncl) { # if the random subset doesn't include all the classes
      return (TRUE) # retry random selection
    }
    
    ## 2. Minimum number of samples per class
    ## Check each class in sub-samples container more than `min_cl` number of samples
    if (any(table(cls) < min_cl)) { # if there is fewer samples
      return (TRUE) # retry random selection
    }
    
    ## 3. Minimum number of unique values per feature
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
    
    ## The `rand_s` satisfies all three quality criteria is returned for LDA
    return (FALSE)
  }


# Perform LDA modeling
# 
# @param data A data frame with z-score values. Rows are samples and columns 
# are features that pass the significance cutoff. 
# @param groups The names of groups for the main class.
# 
ldaFunction <- function (data, groups) {
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
    ss <- data[,-match("class", colnames(data))]
    xy.matrix <- as.matrix(ss) # the original feature matrix

    ## Transform the original feature matrix
    LD <- xy.matrix %*% w.unit # discriminant scores
    
    ## Calculating Effect Size for each feature
    ## Effect Size = the absolute difference between the mean discriminant scores of the two classes
    effect_size <-
        abs(mean(LD[data[, "class"] == 1]) - mean(LD[data[, "class"] == 0]))
    
    ## Coefficient scaling
    ## Scale the unit-normalized LDA coefficients by the effect size
    scal <- w.unit * effect_size # scaled LDA coefficient

    coeff <- vector("numeric", length(scal)) # absolute scaled LDA coefficient
    for (v in seq_along(scal)) {
        if (!is.na(scal[v])) {
            coeff[v] <- abs(scal[v])
        } else{
            coeff[v] <- 0
        }
    }
  
    ## Feature Importance
    lda.means.diff <- (lda.fit$means[2,] - lda.fit$means[1,]) # difference between the class means for each feature
    res <- (lda.means.diff + coeff) / 2
    
    return(res)
}


## Kruskal-Wallis Rank Sum Test for the classes
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
#' between two sample groups for multiple features and uses linear discriminant analysis
#' to establish their effect sizes. Subclass information for each class can be incorporated
#' into the analysis (see examples). Features with large differences between two sample groups
#' are identified as biomarkers.
#'
#' @details
#' The LEfSe method expects relative abundances in the `expr` input. A warning
#' will be emitted if the column sums do not result in 1. Use the \code{relativeAb}
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
#' @param trim.names Default is `FALSE`. If `TRUE`, this function extracts 
#' the most specific taxonomic rank of organism.
#' @param checkAbundances `logical(1)` Whether to check if the assay data in the
#'   `relab` input are relative abundances or counts. If counts are found, a
#'   warning will be emitted (default `TRUE`).
#' @param \ldots Additional inputs to lower level functions (not used).
#' @return
#' The function returns a `data.frame` with two columns, which are
#' names of features and their LDA scores.
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
#'     res_group <- lefser(relativeAb(zeller14), 
#'                         groupCol = "study_condition")
#'     head(res_group)
#'
#'     # (2) Using classes and sublasses
#'     data(zeller14)
#'     # exclude 'adenoma'
#'     zeller14 <- zeller14[, zeller14$study_condition != "adenoma"]
#'     res_block <- lefser(relativeAb(zeller14), 
#'                         groupCol = "study_condition", 
#'                         blockCol = "age_category")
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
      
    ## Check whether relative abundance is provided or not
    if (checkAbundances && !identical(all.equal(colSums(relab_data), 
                                                rep(1e6, ncol(relab_data)), 
                                                check.attributes = FALSE), 
                                      TRUE)) {
        warning("Convert counts to relative abundances with 'relativeAb()'")
    }
     
    ## Extract the class/group information   
    groupf <- colData(relab)[[groupCol]]
    groupf <- as.factor(groupf)
    lgroupf <- levels(groupf)
    if (is.null(groupf) || !identical(length(lgroupf), 2L)) {
        msg <- "'groupCol' must refer to a valid dichotomous (two-level) variable"
        stop(msg) # ensure the class has only two levels
    }
    message(
        "The outcome variable is specified as '", groupCol,
        "' and the reference category is '", lgroupf[1],
        "'.\n See `?factor` or `?relevel` to change the reference category."
    )
    
    ## Kruskal-Wallis Rank Sum Test for the classes
    relab_sub <- filterKruskal(relab = relab_data, 
                               group = groupf, 
                               p.value = kruskal.threshold)

    ## Wilcoxon Rank-Sum Test for the sub-classes
    if (!is.null(blockCol)) {
        block <- as.factor(colData(relab)[[blockCol]])
        block <- droplevels(block)
        ## z-statistics result of the features passing the significance cut-off
        relab_sub <- fillPmatZmat(group = groupf, 
                                  block = block, 
                                  relab_sub = relab_sub, 
                                  p.threshold = wilcox.threshold)
    }
    
    ## Return an empty data table if there is no significant features
    if(nrow(relab_sub) == 0L){
      return(.return_no_results())
    }

    # transposes matrix and add a "class" (i.e., groupf) column
    # matrix converted to dataframe
    relab_sub_t <- t(relab_sub)
    relab_sub_t_df <- as.data.frame(relab_sub_t)
    relab_sub_t_df <- createUniqueValues(df = relab_sub_t_df, group = groupf)
    relab_sub_t_df <- cbind(relab_sub_t_df, class = groupf)

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