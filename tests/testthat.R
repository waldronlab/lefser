library(testthat)
library(lefser)

test_check("lefser")

data(zeller14)

zeller14sub <- zeller14[, zeller14$study_condition != "adenoma"]

zeller14ra <- lefser::relativeAb(zeller14sub)


expect_equal(assay(zeller14ra)$rel_abs[1, 1], assay(zeller14ra)$counts[1, 1] * 1e6)

