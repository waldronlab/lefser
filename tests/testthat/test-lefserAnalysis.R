test_that("lefserAnalysis and lefserPlot work", {
  data("zeller14")
  zellersub <- zeller14[1:200, zeller14$study_condition != "adenoma"]
  expect_error(lefserAnalysis(zellersub),
               "Group assignment 'grp' must be specified")
  # assign '0' class to 'conrol' and '1' to 'CRC' (i.e., colorectal cancer)
  zellersub$GROUP <-
    ifelse(zellersub$study_condition == "control", 0, 1)
  set.seed(1)
  expect_message(results <- lefserAnalysis(zellersub), "Length of block: 0")
  expect_equal(colnames(results), c("Names", "scores"))
  expect_equal(results[1, "Names"], "p__Firmicutes")
  expect_equal(results[nrow(results), "Names"], "o__Bacteroidales")
  expect_equal(results[1, "scores"],-6.431365, tolerance = 1e-4)
  expect_equal(results[nrow(results), "scores"], 6.402655, tolerance = 1e-4)
  # test with blocks: assign '0' class to 'adult' and '1' to 'senior'
  zellersub$BLOCK <- ifelse(zellersub$age_category == "adult", 0, 1)
  expect_message(results2 <- lefserAnalysis(zellersub), "Length of block: 157")
  expect_equal(nrow(results2), 4)
  expect_equal(
    results2$Names,
    c(
      "o__Lactobacillales",
      "s__Ruminococcus_sp_5_1_39BFAA",
      "s__Eubacterium_hallii",
      "s__Streptococcus_salivarius"
    )
  )
  expect_equal(results2$scores,
               c(-5.797806,-5.285158,  3.780228,  4.389330),
               tolerance = 1e-4)
  plt <- lefserPlot(results2)
  expect_s3_class(plt, "ggplot")
})

