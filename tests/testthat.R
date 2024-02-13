library(testthat)
library(lefser)

  

test_check("lefser")

data(zeller14)

zeller14sub <- zeller14[, zeller14$study_condition != "adenoma"]

zeller14ra <- lefser::relativeAb(zeller14sub)


expect_warning(lefser(zeller14sub, groupCol = "study_condition"), regexp = "^Convert counts.+")
expect_equal(assayNames(zeller14ra), c("rel_abs", "exprs"))
expect_equal(colSums(assay(zeller14ra, "rel_abs")), rep(1e6, ncol(zeller14ra)), check.attributes = FALSE)
set.seed(1)
expect_no_warning(x <- lefser(zeller14ra, groupCol = "study_condition"))
