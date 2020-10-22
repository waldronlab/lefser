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
  expect_error(lefser(zellersub, groupCol = NULL, blockCol = NULL))
  results <- withr::with_seed(1,
    lefser(zellersub, groupCol = "study_condition", blockCol = NULL)
  )
  expect_equal(colnames(results), c("Names", "scores"))
  expect_true(checkEnding(results, 1, "Names", "p__Firmicutes`"))
  expect_true(
    checkEnding(results, nrow(results), "Names", "o__Bacteroidales`")
  )
# TODO: compare results between LEfSe and lefser
  expect_equal(results[1, "scores"], -6.431365, tolerance = tol)
  expect_equal(results[nrow(results), "scores"], 6.402655, tolerance = tol)
  results2 <- withr::with_seed(1,
    lefser(
        zellersub, groupCol = "study_condition", blockCol = "age_category"
    )
  )
  expect_equal(colnames(results2), c("Names", "scores"))
  expect_true(checkEnding(results2, 1, "Names", "o__Lactobacillales`"))
  expect_true(checkEnding(
    results2, nrow(results2), "Names", "s__Eubacterium_hallii`"
  ))
  expect_equal(results2[1, "scores"], -5.792313, tolerance = tol)
  expect_equal(results2[nrow(results2), "scores"], 5.20564, tolerance = tol)

  expect_equal(nrow(results2), 3)
  expect_true(all(
    mapply(endsWith, results2$Names, c(
      "o__Lactobacillales`",
      "s__Ruminococcus_sp_5_1_39BFAA`",
      "s__Eubacterium_hallii`"
    ))
  ))
  expect_equal(results2$scores,
               c(-5.792313,-5.079186, 5.205640),
               tolerance = tol)
  plt <- lefserPlot(results2)
  expect_s3_class(plt, "ggplot")
})

