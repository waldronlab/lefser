wilcox_pstats <- function(datamat, class, index) {
  apply(datamat, 1L, function(values) {
    wx <- suppressWarnings({
        coin::wilcox_test(values ~ class, subset = index)
    })
    c(pvalue = coin::pvalue(wx), statistic = coin::statistic(wx))
  })
}

## Wilcoxon Rank-Sum Test for the sub-classes
fillPmatZmat <- function(class,
                         subclass,
                         relab_sub,
                         p.threshold,
                         method)
{
  if (!nrow(relab_sub))
    return(relab_sub)
  # creates a list of boolean vectors, each vector indicates
  # existence (TRUE) or absence (FALSE) of a class/sub-class combination
  combos <- apply(
    expand.grid(levels(class), levels(subclass)), 1L, paste0, collapse = "."
  )
  combined <- interaction(class, subclass)
  whichlist <- lapply(
    setNames(nm = sort(combos)), function(combo) which(combo == combined)
  )

  ## uses Wilcoxon rank-sum test to test for significant differential abundances between
  ## subclasses of one class against subclasses of all other classes; results are saved in
  ## "pval_mat" and "z_mat" matrices
  ssubclass <- seq_along(levels(subclass))
  iters <- expand.grid(ssubclass, ssubclass + length(ssubclass))
  inds <- apply(iters, 1L, function(x) unlist(whichlist[x], use.names = FALSE))

  z_mat <- pval_mat <- matrix(
    NA, nrow = nrow(relab_sub), ncol = length(inds),
    dimnames = list(rownames(relab_sub), NULL)
  )
  for (i in seq_along(inds)) {
    results <- wilcox_pstats(relab_sub, class = class, index = inds[[i]])
    pval_mat[, i] <- stats::p.adjust(results["pvalue", ], method = method)
    z_mat[, i] <- results["statistic", ]
  }

  ## converts "pval_mat" into boolean matrix "logical_pval_mat" where
  ## p-values <= wilcoxon.threshold
  logical_pval_mat <- !is.na(pval_mat) & pval_mat <= p.threshold * 2.0

  ## determines which rows (features) have all p-values<=0.05
  ## and selects such rows from the matrix of z-statistics
  sub <- rowSums(logical_pval_mat) == ncol(logical_pval_mat)
  z_mat_sub <- z_mat[sub, , drop = FALSE]
  ## confirms that z-statistics of a row all have the same sign
  sub <- abs(rowSums(z_mat_sub)) == rowSums(abs(z_mat_sub))
  relab_sub[names(sub[sub]), , drop = FALSE]
}

## ensures that more than half of the values in each for each feature are unique
## if that is not the case then a count value is altered by adding it to a small value
## generated via normal distribution with mean=0 and sd=5% of the count value
createUniqueValues <- function(df, class){
  orderedrows <- rownames(df)
  splitdf <- split(df, class)
  maxim <- vapply(table(class), function(x) max(x * 0.5, 4), numeric(1L))
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


# Perform LDA modeling
#
# @param data A data frame with z-score values. Rows are samples and columns
# are features that pass the significance cutoff.
# @param classes The names of classes for the main class.
#
ldaFunction <- function (data, classes) {

    ## Fitting LDA model
    lda.fit <- lda(class ~ ., data = data)
    w <- lda.fit$scaling[, 1] # extract LDA coefficients
    w.unit <- w / sqrt(sum(w ^ 2)) # scaling LDA coefficient by their Euclidean norm to get unit-normalized coefficient

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
filterKruskal <- function(relab, class, p.value, method = method) {
  # applies "kruskal.test.alt" function to each row (feature) of relab
  # to detect differential abundance between classes, 0 and 1
  kw.res <- apply(relab, 1L, function(x) {
    kruskal.test(x ~ class)[["p.value"]]
  })
  # TRUE for p-values less than or equal to kw.threshold
  kw.res <- stats::p.adjust(kw.res, method = method)
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
#' between two sample classes for multiple features and uses linear discriminant analysis
#' to establish their effect sizes. Subclass information for each class can be incorporated
#' into the analysis (see examples). Features with large differences between two sample classes
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
#' Sum Test (default 0.05). If multiple hypothesis testing is performed, this
#' threshold is applied to corrected p-values.
#' @param wilcox.threshold numeric(1) The p-value for the Wilcoxon Rank-Sum Test
#' when 'subclassCol' is present (default 0.05). If multiple hypothesis testing is
#' performed, this threshold is applied to corrected p-values.
#' @param lda.threshold numeric(1) The effect size threshold (default 2.0).
#' @param classCol character(1) Column name in `colData(relab)` indicating
#' class, usually a factor with two levels (e.g., `c("cases", "controls")`;
#' default "CLASS").
#' @param subclassCol character(1) Optional column name in `colData(relab)`
#' indicating the subclasses, usually a factor with two levels (e.g.,
#' `c("adult", "senior")`; default NULL), but can be more than two levels.
#' @param groupCol (Deprecated) Column name in `colData(relab)` indicating
#'   groups, usually a factor with two levels (e.g., `c("cases", "controls")`;
#'   default "GROUP").
#' @param blockCol (Deprecated) Optional column name in `colData(relab)`
#'   indicating the blocks, usually a factor with two levels (e.g., `c("adult",
#'   "senior")`; default NULL).
#' @param assay The i-th assay matrix in the `SummarizedExperiment` ('relab';
#' default 1).
#' @param trim.names Default is `FALSE`. If `TRUE`, this function extracts
#' the most specific taxonomic rank of organism.
#' @param checkAbundances `logical(1)` Whether to check if the assay data in the
#'   `relab` input are relative abundances or counts. If counts are found, a
#'   warning will be emitted (default `TRUE`).
#' @param method Default is "none" as in the original LEfSe implementation.
#' Character string of length one, passed on to \link[stats]{p.adjust} to set
#' option for multiple testing. For multiple pairwise comparisons, each comparison
#' is adjusted separately. Options are "holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr" (synonym for "BH"), and "none".
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
#' @importFrom testthat capture_warnings
#' @import SummarizedExperiment
#'
#' @examples
#'
#'     data(zeller14)
#'     zeller14 <- zeller14[, zeller14$study_condition != "adenoma"]
#'     tn <- get_terminal_nodes(rownames(zeller14))
#'     zeller14tn <- zeller14[tn,]
#'     zeller14tn_ra <- relativeAb(zeller14tn)
#'
#'     # (1) Using classes only
#'     res_class <- lefser(zeller14tn_ra,
#'                         classCol = "study_condition")
#'     # (2) Using classes and sub-classes
#'     res_subclass <- lefser(zeller14tn_ra,
#'                         classCol = "study_condition",
#'                         subclassCol = "age_category")
#'
#' @export
lefser <-
  function(relab,
           kruskal.threshold = 0.05,
           wilcox.threshold = 0.05,
           lda.threshold = 2.0,
           classCol = "CLASS",
           subclassCol = NULL,
           assay = 1L,
           trim.names = FALSE,
           checkAbundances = TRUE,
           method = "none",
           ...,
           expr,
           groupCol = "GROUP",
           blockCol = NULL
) {
    if (!missing(expr)) {
        .Deprecated(
            msg = "The 'expr' argument is deprecated, use 'relab' instead."
        )
        relab_data <- assay(expr, i = assay)
    } else {
        relab_data <- assay(relab, i = assay)
    }

    if (!missing(groupCol)){
        .Deprecated(
            msg = "The 'groupCol' argument is deprecated, use 'classCol' instead."
        )
        classCol <- groupCol
    }
    if (!missing(blockCol)){
        .Deprecated(
            msg = "The 'blockCol' argument is deprecated, use 'subclassCol' instead."
        )
        subclassCol <- blockCol
    }
    ## Check whether relative abundance is provided or not
    if (checkAbundances && !identical(all.equal(colSums(relab_data),
                                                rep(1e6, ncol(relab_data)),
                                                check.attributes = FALSE),
                                      TRUE)) {
        warning("Convert counts to relative abundances with 'relativeAb()'")
    }

    ## Extract the class/subclass information
    classf <- colData(relab)[[classCol]]
    classf <- as.factor(classf)
    lclassf <- levels(classf)
    if (is.null(classf) || length(lclassf) < 2L)
        stop("'classCol' must be a valid two or more level variable")
    message(
        "The outcome variable is specified as '", classCol,
        "' and the reference category is '", lclassf[1],
        "'.\n See `?factor` or `?relevel` to change the reference category."
    )

    ## Kruskal-Wallis Rank Sum Test for the classes
    relab_sub <- filterKruskal(relab = relab_data,
                               class = classf,
                               p.value = kruskal.threshold,
                               method = method)

    ## Wilcoxon Rank-Sum Test for the sub-classes
    if (!is.null(subclassCol)) {
        subclass <- as.factor(colData(relab)[[subclassCol]])
        subclass <- droplevels(subclass)
        ## z-statistics result of the features passing the significance cut-off
        relab_sub <- fillPmatZmat(class = classf,
                                  subclass = subclass,
                                  relab_sub = relab_sub,
                                  p.threshold = wilcox.threshold,
                                  method = method)
    }

    ## Return an empty data table if there is no significant features
    if(nrow(relab_sub) == 0L){
      return(.return_no_results())
    }

    ## Transposed relative abundance matrix with the 'class' column
    relab_sub_t <- t(relab_sub)
    relab_sub_t_df <- as.data.frame(relab_sub_t)
    # relab_sub_t_df <- createUniqueValues(df = relab_sub_t_df, class = classf)
    relab_sub_t_df <- cbind(relab_sub_t_df, class = classf)

    ## LDA model
    warn <- testthat::capture_warnings(
        raw_lda_scores <- ldaFunction(relab_sub_t_df, lclassf)
    )

    ## Warning collinearity and recommend `get_terminal_nodes`
    if (length(warn) && nzchar(warn)) {
        msg <- "Variables in the input are collinear. Try only with the terminal nodes using `get_terminal_nodes` function"
        warning(msg)
    }

    ## Processing LDA scores
    processed_scores <-
      sign(raw_lda_scores) * log((1 + abs(raw_lda_scores)), 10)
    processed_sorted_scores <- sort(processed_scores)   # sorting of scores
    scores_df <- data.frame(features = names(processed_sorted_scores),
                            scores = as.vector(processed_sorted_scores),
                            stringsAsFactors = FALSE)
    scores_df <- .trunc(scores_df, trim.names)   # short-form of taxa name

    ## Filter with LDA threshold
    threshold_scores <- abs(scores_df$scores) >= lda.threshold
    res_scores <- scores_df[threshold_scores, ]
    class(res_scores) <- c("lefser_df", class(res_scores))
    attr(res_scores, "classes") <- lclassf

    ## Add attributes with argument values
    ## This is used for plottin functions
    attr(res_scores, "inputSE") <- relab
    attr(res_scores, "kth") <-  kruskal.threshold
    attr(res_scores, "wth") <- wilcox.threshold
    attr(res_scores, "ldath") <- lda.threshold
    attr(res_scores, "class") <- classCol
    attr(res_scores, "subclass") <- subclassCol
    attr(res_scores, "method") <- method
    attr(res_scores, "lclassf") <- lclassf[1]
    res_scores
 }

.return_no_results <- function() {
    message("No significant features found.")
    res_scores <- data.frame(features = character(), scores = numeric())
    class(res_scores) <- c("lefser_df", class(res_scores))
    return(res_scores)
}
