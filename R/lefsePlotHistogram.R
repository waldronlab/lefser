.prepareDataHistogram <- function(res) {
    if (!.isLefser(res)) {
        stop(
            "Expected an object of class 'lefser_df'.",
            " Only the output of a lefser function call is valid.",
            call. = FALSE
        )
    }
    tse <- attr(res, "inputSE")
    groupCol <- attr(res, "grp")
    blockCol <- attr(res, "blk")
    selectCols <- c(groupCol, blockCol)
    sampleData <- as.data.frame(SummarizedExperiment::colData(tse))
    sampleData <- sampleData[, selectCols, drop = FALSE]
    sampleData <- tibble::rownames_to_column(sampleData, var = "sample")
    output <- tse[res[["features"]],] |> 
        SummarizedExperiment::assay() |> 
        as.data.frame() |> 
        tibble::rownames_to_column(var = "features") |> 
        tidyr::pivot_longer(
            cols = 2:tidyselect::last_col(),
            names_to = "sample",
            values_to = "abundance"
        ) |> 
        dplyr::left_join(sampleData, by = "sample")
    return(output)
}

.isLefser <- function(x) {
    "lefser_df" %in% class(x)
}