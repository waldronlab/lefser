#' Plot Feature
#'
#' \code{lefserPlotFeat} plots the abundance data of a DA feature across all
#' samples.
#'
#' @param res An object of class lefser_df,
#' output of the \code{lefser} function.
#' @param fName A character string. The name of a feature in the lefser_df
#' object.
#' @param colors Colors corresponding to class 0 and 1.
#' Options: "c" (colorblind), "l" (lefse), "g" (greyscale). Defaults to "c".
#' This argument also accepts a character(2) with two color names.
#'
#' @details
#' The solid lines represent the mean by class or by class+subclass
#' (if the subclass variable is present).
#' The dashed lines represent the median by class or by class+subclass
#' (if the subclass variable is present).
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#'
#' data(zeller14)
#' zeller14 <- zeller14[, zeller14$study_condition != "adenoma"]
#' tn <- get_terminal_nodes(rownames(zeller14))
#' zeller14tn <- zeller14[tn,]
#' zeller14tn_ra <- relativeAb(zeller14tn)
#'
#' # (1) Using classes only
#' res_class <- lefser(zeller14tn_ra,
#'                     classCol = "study_condition")
#' # (2) Using classes and sub-classes
#' res_subclass <- lefser(zeller14tn_ra,
#'                     classCol = "study_condition",
#'                     subclassCol = "age_category")
#' plot_class <- lefserPlotFeat(res_class, res_class$features[[1]])
#' plot_subclass <- lefserPlotFeat(res_subclass, res_subclass$features[[2]])
#'
lefserPlotFeat <- function(res, fName, colors = "colorblind") {
    dat <- .prepareDataHistogram(res = res, fName = fName)
    refclass <- attr(res, "lclassf")
    vLinePos <- which(dat$classCol != refclass)[1] - 0.5
    colVar <- .selectPalette(colors)
    l <- split(dat, dat$classCol)
    maxYVal <- max(dat$abundance)
    vals <- purrr::map(l, ~ {
        head(.x[["sample"]], 1)
    })
    sumDat <- .sumDatHist(dat)
    cond <- "subclassCol" %in% colnames(dat)
    if (isFALSE(cond)) {
        p <- dat |>
            ggplot2::ggplot(
                data = dat, mapping = ggplot2::aes(sample, abundance)
            ) +
            ggplot2::geom_col(
                mapping = ggplot2::aes(fill = classCol), width = 1
            ) +
            ggplot2::scale_fill_manual(
                values = colVar
            )
    } else if (isTRUE(cond)) {
        p <- dat |>
            ggplot2::ggplot(
                data = dat, mapping = ggplot2::aes(sample, abundance)
            ) +
            ggplot2::geom_col(
                mapping = ggplot2::aes(fill = subclassCol), width = 1
            ) +
            ggplot2::scale_fill_manual(
                values = colVar
            )
    }
    p <- p +
        ggplot2::scale_y_continuous(
            labels = scales::scientific, expand = c(0, 0)
        ) +
        ggplot2::labs(
            x = "Samples", y = "Relative abundance",
            title = stringr::str_wrap(fName, whitespace_only = FALSE)
        ) +
        ggplot2::annotate(
            geom = "label",
            x = vals[[levels(dat$classCol)[1]]],
            y = maxYVal,
            label = levels(dat$classCol)[1],
            hjust = 0, vjust = 1
        ) +
        ggplot2::annotate(
            geom = "label",
            x = vals[[levels(dat$classCol)[2]]],
            y = maxYVal,
            label = levels(dat$classCol)[2],
            hjust = 0, vjust = 1
        ) +
        ggplot2::geom_vline(xintercept = vLinePos) +
        theme_bw() +
        theme(
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_blank(),
            legend.title = ggplot2::element_blank()
        )
    for (i in seq_along(sumDat)) {
        p <- p +
            ggplot2::geom_segment(
                data = sumDat[[i]],
                mapping = ggplot2::aes(
                    x = x1, xend = x2,
                    y = mean1, yend = mean2
                ),
                linetype = 1
            ) +
            ggplot2::geom_segment(
                data = sumDat[[i]],
                mapping = ggplot2::aes(
                    x = x1, xend = x2,
                    y = median1, yend = median2
                ),
                linetype = 2
            )
    }
    return(p)
}

.sumDatHist <- function(x) {
    cond <- "subclassCol" %in% colnames(x)
    if (isFALSE(cond)) {
        l <- split(x, x$classCol)
    } else if (isTRUE(cond)) {
       l <- x |>
           dplyr::mutate(
               classes = paste0(.data$classCol, "--", .data$subclassCol)
           ) |>
           {\(y) split(y, y[["classes"]])}()
    }
    sumDat <- purrr::map(l, ~ {
        meanVal <- mean(.x[["abundance"]] )
        medianVal <- median(.x[["abundance"]] )
        minVal <- head(.x[["sample"]], 1)
        maxVal <- tail(.x[["sample"]], 1)
        data.frame(
            x1 = minVal, x2 = maxVal, mean1 = meanVal, mean2 = meanVal,
            median1 = medianVal, median2 = medianVal
        )
    })
    return(sumDat)
}

.prepareDataHistogram <- function(res, fName) {
    if (!.isLefser(res)) {
        stop(
            "Expected an object of class 'lefser_df'.",
            " Only the output of a lefser function call is valid.",
            call. = FALSE
        )
    }
    tse <- attr(res, "inputSE")
    classCol <- attr(res, "class_arg")
    subclassCol <- attr(res, "subclass_arg")
    refclass <- attr(res, "lclassf")
    selectCols <- c(classCol, subclassCol)
    sampleData <- as.data.frame(SummarizedExperiment::colData(tse))
    sampleData <- sampleData[, selectCols, drop = FALSE]
    sampleData <- tibble::rownames_to_column(sampleData, var = "sample")
    dat <- tse[res[["features"]],] |>
        SummarizedExperiment::assay() |>
        as.data.frame() |>
        tibble::rownames_to_column(var = "features") |>
        tidyr::pivot_longer(
            cols = 2:tidyselect::last_col(),
            names_to = "sample",
            values_to = "abundance"
        ) |>
        dplyr::left_join(sampleData, by = "sample")
    dat[[classCol]] <- forcats::fct_relevel(dat[[classCol]], refclass)

    if (is.null(subclassCol)) {
        dat <- dat |>
            dplyr::filter(.data[["features"]] == fName) |>
            dplyr::arrange(
                .data[[classCol]],
                dplyr::desc(.data[["abundance"]])
            ) |>
            dplyr::mutate(sample = forcats::fct_inorder(.data[["sample"]]))
        colnames(dat)[which(colnames(dat) == classCol)] <- "classCol"
    } else if (!is.null(subclassCol)) {
        dat <- dat |>
            dplyr::filter(.data[["features"]] == fName) |>
            dplyr::arrange(
                .data[[classCol]],
                dplyr::desc(.data[[subclassCol]]),
                dplyr::desc(.data[["abundance"]])
            ) |>
            dplyr::mutate(sample = forcats::fct_inorder(.data[["sample"]]))
        colnames(dat)[which(colnames(dat) == classCol)] <- "classCol"
        colnames(dat)[which(colnames(dat) == subclassCol)] <- "subclassCol"
    }
    return(dat)
}

.isLefser <- function(x) {
    "lefser_df" %in% class(x)
}
