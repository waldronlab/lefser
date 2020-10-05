createGroupBlockMatrixGroups <-
  function(se, groupCol, blockCol, assay) {
    expr <- assay(se, i = assay)
    grp <- colData(se)[[groupCol]]
    blk <- colData(se)[[blockCol]]
    
    if (is.null(grp))
      stop("Group assignment 'grp' must be specified")
    grp <- as.factor(grp)
    groups <- levels(grp)
    if (length(groups) != 2L)
      stop(
        "Group classification is not dichotomous:\n",
        "Found (",
        paste(groups, collapse = ", "),
        ")"
      )
    list(group = grp,
         block = as.factor(blk),
         expr = expr)
  }

fillPmatZmat <- function(group,
                         block,
                         expr_sub,
                         wilcoxon.threshold)
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
      logical.list.of.subclasses <- 
        append(logical.list.of.subclasses,
               list(group == sort(unique(group))[i] &
                      block == sort(unique(block))[j]))
    }
  }
  
  ## creates empty matrices that will be filled by p-values and z-statistics from 
  ## Wilcoxon rank-sum test
  pval_mat <-
    matrix(0, nrow(expr_sub), length(unique(block)) ^ length(unique(group)))
  z_mat <-
    matrix(0, nrow(expr_sub), length(unique(block)) ^ length(unique(group)))
  
  rownames(pval_mat) <- rownames(expr_sub)
  rownames(z_mat) <- rownames(expr_sub)
  
  ## uses Wilcoxon rank-sum test to test for significant differential abundances between
  ## subclasses of one class against subclasses of all othe classes; results are saved in
  ## "pval_mat" and "z_mat" matrices
  c <- 0
  for (i in seq_along(unique(block))) {
    for (j in seq_along(unique(block))) {
      j <- j + length(unique(block))
      mat = cbind(expr_sub[, logical.list.of.subclasses[[i]]],
                  expr_sub[, logical.list.of.subclasses[[j]]])
      group_for_mat = factor(c(rep(
        0, sum(logical.list.of.subclasses[[i]])
      ),
      rep(
        1, sum(logical.list.of.subclasses[[j]])
      )))
      pval_mat[, c + 1] <-
        suppressWarnings(apply(mat, 1, wilcox_test_pvalue, group = group_for_mat))
      z_mat[, c + 1] <-
        suppressWarnings(apply(mat, 1, wilcox_test_z, group = group_for_mat))
      c = c + 1
    }
  }
  
  ## number of p-values per feature
  num_comp <- length(unique(block)) ^ length(unique(group))
  
  ## converts "pval_mat" into boolean matrix "logical_pval_mat" where 
  ## p-values <= wilcoxon.threshold
  logical_pval_mat <- pval_mat <= wilcoxon.threshold
  
  ## determines which rows (features) have all p-values<=0.05
  ## and selects such rows from the matrix of z-statistics
  sub <- which(rowSums(logical_pval_mat) == num_comp)
  z_mat_sub <- z_mat[sub,]
  
  # confirms that z-statistics of a row all have the same sign
  sub <- abs(rowSums(z_mat_sub)) == abs(rowSums(abs(z_mat_sub)))
  expr_sub <- expr_sub[names(sub[sub]), ]
}

## ensures that more than half of the values in each for each feature are unique
## if that is not the case then a count value is altered by adding it to a small value
## generated via normal distribution with mean=0 and sd=5% of the count value
createUniqueValues <- function(expr_sub_t_df, groups) {
  df <- expr_sub_t_df
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
  expr_sub_t_df <- df
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
#' groups, usually a factor with two levels (e.g., `c("cases", "controls")`).
#' @param blockCol character(1) Column name in `colData(expr)` indicating the
#' blocks, usually a factor with two levels (e.g., `c("adult", "senior")`).
#' @param assay The i-th assay matrix in the `SummarizedExperiment` ('expr').
#' Defaults to 1.
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
#' @export
lefserAnalysis <-
  function(expr,
           kw.threshold = 0.05,
           wilcoxon.threshold = 0.05,
           lda.threshold = 2.0,
           groupCol = "GROUP",
           blockCol = "BLOCK",
           assay = 1L)
  {
    groupBlockMatrixGroups <-
      createGroupBlockMatrixGroups(expr, groupCol, blockCol, assay)
    group <- groupBlockMatrixGroups$group
    block <- groupBlockMatrixGroups$block
    expr <- groupBlockMatrixGroups$expr
    groups <- levels(group)
    
    # extracts p-values from Kruskal-Wallis Rank Sum Test
    kruskal.test.alt <- function(x, group) {
      kruskal.test(x ~ group)$p.value
    }
    
    # applies "kruskal.test.alt" function to each row (feature) of expr
    # to detect differential abundance between classes, 0 and 1
    kw.res <- apply(expr, 1, kruskal.test.alt, group = group)
    
    # selects p-values less than or equal to kw.threshold
    kw.sub <- kw.res <= kw.threshold
    
    # eliminates NAs
    kw.sub[is.na(kw.sub)] <- FALSE
    
    # extracts features with statistically significant differential abundance
    # from "expr" matrix
    expr_sub <- expr[kw.sub,]
    message("Length of block: ", length(block))
    if (length(block) != 0) {
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
    min_cl <-
      as.integer(min(table(expr_sub_t_df$class)) * 2 / 3 * 2 / 3 *
                   0.5)
    # if min_cl is less than 1, then make it equal to 1
    min_cl <- max(min_cl, 1)
    
    # lda_fn repeated 30 times, producing a matrix of 30 scores per feature
    eff_size_mat <-
      replicate(30, suppressWarnings(ldaFunction(
        expr_sub_t_df, lfk, rfk, min_cl, ncl, groups
      )), simplify = T)
    
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
