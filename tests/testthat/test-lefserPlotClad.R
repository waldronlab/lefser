data("zeller14")
z14 <- zeller14[, zeller14$study_condition != "adenoma"]
tn <- get_terminal_nodes(rownames(z14))
z14tn <- z14[tn, ]
z14tn_ra <- relativeAb(z14tn)

res <- lefser(z14tn_ra, groupCol = "study_condition")
resAll <- lefserAllRanks(relab = z14tn_ra, groupCol = "study_condition")

test_that("lefserAllRanks works", {
    expect_s3_class(resAll, "lefser_df_all")
})

test_that("lefserPlotClad works", {
    ggt <- lefserPlotClad(df = resAll)
    expect_s3_class(ggt, "ggtree")
})