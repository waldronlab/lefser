data("zeller14")
z14 <- zeller14[, zeller14$study_condition != "adenoma"]
tn <- get_terminal_nodes(rownames(z14))
z14tn <- z14[tn, ]
z14tn_ra <- relativeAb(z14tn)

res <- lefser(z14tn_ra, classCol = "study_condition")

z14_input <- rowNames2RowData(z14tn_ra)
test_that("lefserClades works", {
    resClades <<- expect_no_warning(
        withCallingHandlers(
            lefserClades(z14_input, classCol = "study_condition"),
            warning = function(w) {
                # TODO: Remove this handler once `mia` package updates `sumCountsAcrossFeatures` 
                # to `aggregateAcrossFeatures` internally and the deprecation warning is resolved.
                if (grepl("sumCountsAcrossFeatures|aggregateAcrossFeatures", w$message)) {
                    invokeRestart("muffleWarning")
                }
            }
        )
    )
    expect_s3_class(resClades, "lefser_df_clades")
})

test_that("lefserPlotClad works", {
    ggt <- expect_no_warning(
        lefserPlotClad(df = resClades)
    )
    expect_s3_class(ggt, "ggtree")
})