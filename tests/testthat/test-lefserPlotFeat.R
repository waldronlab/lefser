
data(zeller14)
zeller14 <- zeller14[, zeller14$study_condition != "adenoma"]
tn <- get_terminal_nodes(rownames(zeller14))
zeller14tn <- zeller14[tn,]
zeller14tn_ra <- relativeAb(zeller14tn)

# (1) Using classes only
res_class <- lefser(zeller14tn_ra,
                    classCol = "study_condition")
# (2) Using classes and sub-classes
res_subclass <- lefser(zeller14tn_ra,
                    classCol = "study_condition",
                    subclassCol = "age_category")
plot_class <- lefserPlotFeat(res_class, res_class$features[[1]])
plot_subclass <- lefserPlotFeat(res_subclass, res_subclass$features[[2]])

test_that("lefserPlotFeat works", {
  expect_s3_class(plot_class, "ggplot")
  expect_s3_class(plot_subclass, "ggplot")
})
