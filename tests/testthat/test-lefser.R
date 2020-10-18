test_that("lefser and lefserPlot work", {
  data("zeller14")
  zellersub <- zeller14[1:150, zeller14$study_condition != "adenoma"]
  expect_error(lefser(zellersub, groupCol = NULL, blockCol = NULL),
               "Group assignment 'groupCol' must be specified")
  set.seed(1)
  results <- lefser(zellersub, groupCol = "study_condition", blockCol = NULL)
  expect_equal(colnames(results), c("Names", "scores"))
  expect_equal(results[1, "Names"], "p__Firmicutes")
  expect_equal(results[nrow(results), "Names"], "o__Bacteroidales")
  expect_equal(results[1, "scores"],-6.431365, tolerance = 0.2)
  expect_equal(results[nrow(results), "scores"], 6.402655, tolerance = 0.2)
  results2 <- lefser(zellersub, groupCol = "study_condition", blockCol = "age_category")
  
  expect_equal(colnames(results2), c("Names", "scores"))
  expect_equal(results2[1, "Names"], "o__Lactobacillales")
  expect_equal(results2[nrow(results2), "Names"], "s__Eubacterium_hallii")
  expect_equal(results2[1, "scores"],-5.792313, tolerance = 0.2)
  expect_equal(results2[nrow(results2), "scores"], 5.20564, tolerance = 0.2)

  expect_equal(nrow(results2), 3)
  expect_equal(
    results2$Names,
    c(
      "o__Lactobacillales",
      "s__Ruminococcus_sp_5_1_39BFAA",
      "s__Eubacterium_hallii"
    )
  )
  expect_equal(results2$scores,
               c(-5.792313,-5.079186, 5.205640),
               tolerance = 0.2)
  plt <- lefserPlot(results2)
  expect_s3_class(plt, "ggplot")
})

