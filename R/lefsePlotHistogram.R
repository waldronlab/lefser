lefsePlotHistogram <- function(
        res, tse, groupCol = "GROUP", blockCol = NULL
) {
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
