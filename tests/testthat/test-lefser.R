checkEnding <- function(resObj, index, column, value) {
    endsWith(resObj[index, column], value)
}

test_that("lefser and lefserPlot work", {
  skip_if_not(identical(Sys.getlocale("LC_COLLATE"), "en_US.UTF-8"))
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
  expect_true(checkEnding(results2, 1, "Names", "g__Ruminococcus`"))
  expect_true(checkEnding(
    results2, nrow(results2), "Names", "f__Clostridiales_noname`"
  ))
  expect_equal(results2[1, "scores"], -5.792313, tolerance = tol)
  expect_equal(results2[nrow(results2), "scores"], 5.20564, tolerance = tol)

  expect_equal(nrow(results2), 6)
  expect_true(all(
    mapply(endsWith, results2$Names, c(
      "g__Ruminococcus`",
      "o__Lactobacillales`",
      "c__Bacilli`",
      "s__Ruminococcus_sp_5_1_39BFAA`",
      "s__Eubacterium_hallii`",
      "f__Clostridiales_noname`"
    ))
  ))

  expect_equal(results2$scores,
    c(-6.1384025872238, -5.90688120617861, -5.89755502841847,
    -5.57735644006024, -5.42866570267924, 4.64159299932332),
    tolerance = tol)
  plt <- lefserPlot(results2)
  expect_s3_class(plt, "ggplot")
})

