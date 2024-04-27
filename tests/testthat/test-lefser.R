checkEnding <- function(resObj, index, column, value) {
    endsWith(resObj[index, column], value)
}

test_that("lefser and lefserPlot work", {
    dataenv <- new.env(parent = emptyenv())
    data("zeller14", package = "lefser", envir = dataenv)
    zeller14 <- dataenv[["zeller14"]]
    
    tol <- 0.2
    zellersub <- zeller14[1:150, zeller14$study_condition != "adenoma"]
    ## subsetting a DataFrame with NULL
    expect_warning(expect_error(lefser(zellersub, groupCol = NULL, blockCol = NULL)))
    zellersub <- relativeAb(zellersub)
    results <- withr::with_seed(1,
                                lefser(zellersub, groupCol = "study_condition", blockCol = NULL)
    )
    # TODO: compare results between LEfSe and lefser
    expect_equal(results[1, "scores"], -3.688676, tolerance = tol)
    expect_equal(results[nrow(results), "scores"], 3.815295, tolerance = tol)
    results2 <- withr::with_seed(1,
                                 lefser(
                                     zellersub, groupCol = "study_condition", blockCol = "age_category"
                                 )
    )
    expect_equal(results2[1, "scores"], -3.58, tolerance = tol)
    expect_equal(results2[nrow(results2), "scores"], 3.83, tolerance = tol)
    expect_equal(nrow(results2), 14)
    
    # Perform text-based checks only if system locale is en_US.UTF-8
    skip_if_not(identical(Sys.getlocale("LC_COLLATE"), "en_US.UTF-8"))
    expect_equal(results2$scores,
                 c(-3.58274695616907, -3.32856514474989, -3.31328636163603, -3.28470366735756, 
                   -2.99971715869742, -2.82238742034892, 2.36885229314682, 2.36885229314682, 
                   2.36885229314682, 2.53194431507908, 2.53194431507908, 3.83089005744698, 
                   3.83093344936543, 3.83093344936543),
                 tolerance = tol)
    
    expect_equal(colnames(results), c("Names", "scores"))
    expect_true(checkEnding(results, 1, "Names", "p__Firmicutes"))
    expect_true(
        checkEnding(results, nrow(results), "Names", "o__Bacteroidales")
    )
    expect_equal(colnames(results2), c("Names", "scores"))
    expect_true(checkEnding(results2, 1, "Names", "g__Ruminococcus"))
    expect_true(checkEnding(
        results2, nrow(results2), "Names", "o__Bacteroidales"
    ))
    expect_true(all(
        mapply(endsWith, results2$Names[1:7], c(
            "g__Ruminococcus",
            "o__Lactobacillales",
            "c__Bacilli",
            "f__Streptococcaceae",
            "s__Ruminococcus_sp_5_1_39BFAA",
            "s__Eubacterium_hallii",
            "c__Deltaproteobacteria"
        ))
    ))
    
    plt <- lefserPlot(results2)
    expect_s3_class(plt, "ggplot")
})

test_that("no significant results behaviors are consistent", {
    suppressPackageStartupMessages(library(lefser))
    dataenv <- new.env(parent = emptyenv())
    data("zeller14", package = "lefser", envir = dataenv)
    zeller14 <- dataenv[["zeller14"]]
    zellersub <- zeller14[1:150, zeller14$study_condition != "adenoma"]
    zellersub <- relativeAb(zellersub)
    expect_message(
        res1 <- lefser(zellersub, groupCol = "study_condition", blockCol = "age_category", wilcox.threshold = 1e-10),
        "No significant features found."
    )
    expect_message(
        res2 <- lefser(zellersub, groupCol = "study_condition", blockCol = "age_category", kruskal.threshold = 1e-10),
        "No significant features found."
    )
    expect_message(
        res3 <- lefser(zellersub, groupCol = "study_condition", blockCol = "age_category", lda.threshold = 1e6),
        "No significant features found."
    )
    expect_equal(nrow(res1), 0L)
    expect_equal(nrow(res2), 0L)
    expect_equal(nrow(res3), 0L)
    expect_equal(colnames(res1), c("Names", "scores"))
    expect_equal(colnames(res2), c("Names", "scores"))
    expect_equal(colnames(res3), c("Names", "scores"))
})

test_that("Relative abundance", {
    data(zeller14)
    zeller14sub <- zeller14[, zeller14$study_condition != "adenoma"]
    zeller14ra <- lefser::relativeAb(zeller14sub)
    
    expect_warning(lefser(zeller14sub, groupCol = "study_condition"), regexp = "^Convert counts.+")
    expect_equal(assayNames(zeller14ra), c("rel_abs", "exprs"))
    expect_equal(colSums(assay(zeller14ra, "rel_abs")), rep(1e6, ncol(zeller14ra)), check.attributes = FALSE)
    set.seed(1)
    expect_warning(x <- lefser(zeller14ra, groupCol = "study_condition"),
                   "variables are collinear")
})