checkEnding <- function(resObj, index, column, value) {
    endsWith(resObj[index, column], value)
}

library(lefser)
dataenv <- new.env(parent = emptyenv())
data("zeller14", package = "lefser", envir = dataenv)
zeller14 <- dataenv[["zeller14"]]
zeller142 <- zeller14[, zeller14$study_condition != "adenoma"]
tn <- get_terminal_nodes(rownames(zeller142))
zellersub <- zeller142[tn, ]
zellersubra <- relativeAb(zellersub)

## Set the level for class
colData(zellersubra)$study_condition <- factor(colData(zellersubra)$study_condition,
                                               levels = c("control", "CRC"))

test_that("lefser and lefserPlot work", {
    
    ## This is a fractional tolerance, i.e. 0.001 -> 0.1% fractional difference allowed
    tol <- 0.001 
    
    ## Subsetting a DataFrame with NULL
    set.seed(1234)
    expect_warning(expect_error(lefser(zellersub, classCol = NULL, subclassCol = NULL)))
    withr::with_seed(1,
      results <- lefser(zellersubra, classCol = "study_condition", subclassCol = NULL)
    )
    # TODO: compare results between LEfSe and lefser
    expect_equal(results[1, "scores"], -4.313678, tolerance = tol)
    expect_equal(results[nrow(results), "scores"], 4.099828, tolerance = tol)
    withr::with_seed(1, 
      results2 <- lefser(zellersubra, classCol = "study_condition", subclassCol = "age_category")
    )
    expect_equal(nrow(results2), 15)
    expect_equal(results2$scores,
                 c(-3.795320, -3.700813, -3.521247, -3.154544, -3.126245,
                   -3.082713, -2.962662,  2.141700,  2.286564,  2.431915,
                   2.579108,  2.706471, 2.834329,  2.941230,  3.336170), 
                 tolerance = tol)
    
    # Perform text-based checks only if system locale is en_US.UTF-8
    skip_if_not(identical(Sys.getlocale("LC_COLLATE"), "en_US.UTF-8"))
    expect_equal(colnames(results), c("features", "scores"))
    expect_true(checkEnding(results, 1, "features", "t__GCF_000209875"))
    expect_true(
        checkEnding(results, nrow(results), "features", "s__Subdoligranulum_unclassified")
    )
    expect_equal(colnames(results2), c("features", "scores"))
    expect_true(checkEnding(results2, 1, "features", "t__GCF_000159975"))
    expect_true(checkEnding(
        results2, nrow(results2), "features", "s__Oscillibacter_unclassified"
    ))
    expect_true(all(
        mapply(endsWith, results2$features[14:15], c(
            "t__GCF_000147675",
            "s__Oscillibacter_unclassified"
        ))
    ))
    
    plt <- lefserPlot(results2)
    expect_s3_class(plt, "ggplot")
    
    expect_error(
        lefserPlot(results2, label.font.face = "invalid"),
        "'arg' should be one of"
    )
})

test_that("no significant results behaviors are consistent", {
    expect_message(
        res1 <- lefser(zellersubra, classCol = "study_condition", subclassCol = "age_category", wilcox.threshold = 1e-10),
        "No significant features found."
    )
    expect_message(
        res2 <- lefser(zellersubra, classCol = "study_condition", subclassCol = "age_category", kruskal.threshold = 1e-10),
        "No significant features found."
    )
    expect_message(
        res3 <- lefser(zellersubra, classCol = "study_condition", subclassCol = "age_category", lda.threshold = 1e6),
        "No significant features found."
    )
    expect_equal(nrow(res1), 0L)
    expect_equal(nrow(res2), 0L)
    expect_equal(nrow(res3), 0L)
    expect_equal(colnames(res1), c("features", "scores"))
    expect_equal(colnames(res2), c("features", "scores"))
    expect_equal(colnames(res3), c("features", "scores"))
})

test_that("Relative abundance", {
    expect_warning(lefser(zellersub, classCol = "study_condition"), regexp = "^Convert counts.+")
    expect_equal(assayNames(zellersubra), c("rel_abs", "exprs"))
    expect_equal(colSums(assay(zellersubra, "rel_abs")), rep(1e6, ncol(zellersubra)), check.attributes = FALSE)
    set.seed(1)
    expect_warning(x <- lefser(zeller142, classCol = "study_condition"),
                   "Variables in the input are collinear")
})

test_that("ldaFunction correctly identifies classes and calculates scores", {
    class_A_data <- data.frame(
        feature1 = c(10, 11, 12),
        feature2 = c(1, 2, 3),
        class = "A"
    )

    class_B_data <- data.frame(
        feature1 = c(1, 2, 3),
        feature2 = c(10, 11, 12),
        class = "B"
    )

    test_data <- rbind(class_A_data, class_B_data)
    test_data$class <- factor(test_data$class, levels = c("A", "B"))

    classes_levels <- levels(test_data$class)

    ## createUniqueValues adds jitter when needed (required for collinear test data);
    ## use a seed for reproducibility; suppressWarnings() hides the LDA
    ## "variables are collinear" warning that may appear for low-sample data
    withr::with_seed(1, {
        lda_scores <- lefser:::ldaFunction(
            data = test_data, classes = classes_levels
        ) |> suppressWarnings()
    })

    expect_type(lda_scores, "double")
    expect_named(lda_scores, c("feature1", "feature2"))
    ## feature1 is enriched in class A (reference), so its score is negative
    ## feature2 is enriched in class B, so its score is positive
    expect_true(lda_scores[["feature1"]] < 0)
    expect_true(lda_scores[["feature2"]] > 0)
})
