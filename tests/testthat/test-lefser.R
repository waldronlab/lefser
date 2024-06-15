checkEnding <- function(resObj, index, column, value) {
    endsWith(resObj[index, column], value)
}

dataenv <- new.env(parent = emptyenv())
data("zeller14", package = "lefser", envir = dataenv)
zeller14 <- dataenv[["zeller14"]]
zeller142 <- zeller14[, zeller14$study_condition != "adenoma"]
tn <- get_terminal_nodes(rownames(zeller142))
zellersub <- zeller142[tn, ]
zellersubra <- relativeAb(zellersub)

test_that("lefser and lefserPlot work", {
    tol <- 0.001 # this is a fractional tolerance, ie 0.001 -> 0.1% fractional difference allowed
    ## subsetting a DataFrame with NULL
    expect_warning(expect_error(lefser(zellersub, groupCol = NULL, blockCol = NULL)))
    withr::with_seed(1,
      results <- lefser(zellersubra, groupCol = "study_condition", blockCol = NULL)
    )
    # TODO: compare results between LEfSe and lefser
    expect_equal(results[1, "scores"], -4.313678, tolerance = tol)
    expect_equal(results[nrow(results), "scores"], 4.099828, tolerance = tol)
    withr::with_seed(1, 
      results2 <- lefser(zellersubra, groupCol = "study_condition", blockCol = "age_category")
    )
    expect_equal(nrow(results2), 15)
        expect_equal(results2$scores,
                 c(-3.79532040165821, -3.70081304070088, -3.52124711962994, -3.35646736354917, 
                   -3.12624477369647, -3.08271306137103, -3.07415263925044, 2.14169975718221, 
                   2.48592299433808, 2.63986709426167, 2.76738470501209, 2.9270111481154, 
                   3.09767990011079, 3.23928501499009, 3.33616975287112)
                 ,
                 tolerance = tol)
    
    # Perform text-based checks only if system locale is en_US.UTF-8
    skip_if_not(identical(Sys.getlocale("LC_COLLATE"), "en_US.UTF-8"))
    expect_equal(colnames(results), c("Names", "scores"))
    expect_true(checkEnding(results, 1, "Names", "t__GCF_000209875"))
    expect_true(
        checkEnding(results, nrow(results), "Names", "s__Subdoligranulum_unclassified")
    )
    expect_equal(colnames(results2), c("Names", "scores"))
    expect_true(checkEnding(results2, 1, "Names", "t__GCF_000159975"))
    expect_true(checkEnding(
        results2, nrow(results2), "Names", "s__Oscillibacter_unclassified"
    ))
    expect_true(all(
        mapply(endsWith, results2$Names[14:15], c(
            "t__GCF_000147675",
            "s__Oscillibacter_unclassified"
        ))
    ))
    
    plt <- lefserPlot(results2)
    expect_s3_class(plt, "ggplot")
})

test_that("no significant results behaviors are consistent", {
    expect_message(
        res1 <- lefser(zellersubra, groupCol = "study_condition", blockCol = "age_category", wilcox.threshold = 1e-10),
        "No significant features found."
    )
    expect_message(
        res2 <- lefser(zellersubra, groupCol = "study_condition", blockCol = "age_category", kruskal.threshold = 1e-10),
        "No significant features found."
    )
    expect_message(
        res3 <- lefser(zellersubra, groupCol = "study_condition", blockCol = "age_category", lda.threshold = 1e6),
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
    expect_warning(lefser(zellersub, groupCol = "study_condition"), regexp = "^Convert counts.+")
    expect_equal(assayNames(zellersubra), c("rel_abs", "exprs"))
    expect_equal(colSums(assay(zellersubra, "rel_abs")), rep(1e6, ncol(zellersubra)), check.attributes = FALSE)
    set.seed(1)
    expect_warning(x <- lefser(zeller142, groupCol = "study_condition"),
                   "variables are collinear")
})

