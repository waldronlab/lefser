
#' Run lefser at different clades
#' 
#' \code{lefesrCaldes} Agglomerates the features abundance at different
#' taxonomic ranks using \code{mia::splitByRanks} and performs lefser at
#' each rank. The analysis is run at the species, genus, family, order,
#' class, and phylum levels.
#'
#' @param relab A (Tree) SummarizedExperiment with full taxonomy in the rowData.
#' @param ... Arguments passed to the \code{lefser} function.
#'
#' @return An object of class 'lefser_df_clades', "lefser_df", and 'data.frame'.
#' 
#' @details
#' 
#' When running \code{lefserClades}, all features with NAs in the rowData will
#' be dropped. This is to avoid creating artificial clades with NAs.
#' 
#' @export
#'
#' @examples
#' 
#' data("zeller14")
#' z14 <- zeller14[, zeller14$study_condition != "adenoma"]
#' tn <- get_terminal_nodes(rownames(z14))
#' z14tn <- z14[tn, ]
#' z14tn_ra <- relativeAb(z14tn)
#' z14_input <- rowNames2RowData(z14tn_ra)
#'
#' resCl <- lefserClades(relab = z14_input, classCol = "study_condition")
#' 
lefserClades <- function(relab, ...) {
    se <- .selectTaxRanks(relab)
    se <- .appendRankLetter(se)
    l <- .dropFeatures(se)
    se <- l[["se"]]
    pathStrings <- l[["pathStrings"]]
    seL <- as.list(mia::splitByRanks(se))
    ## Kingdom would not be informative
    # seL <- seL[!names(seL) %in% "kingdom"]
    msgRanks <- paste(names(seL), collapse = ", ")
    msgRanks <- msgRanks[length(msgRanks):1]
    message(
        "lefser will be run at the ", msgRanks, " level."
    )
    seL <- purrr::imap(seL, function(x, idx) {
        seVar <- x
        row_data <- as.data.frame(SummarizedExperiment::rowData(seVar))
        row_data <- row_data[,1:which(colnames(row_data) == idx), drop = FALSE]
        # row_data <- purrr::discard(row_data, ~ all(is.na(.x)))
        SummarizedExperiment::rowData(seVar) <- S4Vectors::DataFrame(row_data)
        newRowNames <- .rowData2PathStrings(seVar)
        BiocGenerics::rownames(seVar) <- newRowNames
        seVar
    })
    resL <- purrr::imap(seL, function(x, idx, ...) {
        message(
            "\n>>>> Running lefser at the ", idx, " level.", 
            " <<<<"
        )
        withCallingHandlers(
            lefser(relab = x,...),
            warning = function(w) {
                if (grepl("relativeAb", w$message)) {
                    invokeRestart("muffleWarning")
                }
            }
        )
    },
    ...
    )
    controlVar <- resL |> 
        purrr::map(~ attr(.x, "lclassf")) |> 
        unlist(use.names = FALSE) |> 
        unique()
    caseVar <- resL |> 
        purrr::map(~ attr(.x, "case")) |> 
        unlist(use.names = FALSE) |> 
        unique()
    classArgVar <- resL |> 
        purrr::map(~ attr(.x, "class_arg")) |> 
        unlist(use.names = FALSE) |> 
        unique()
    subclassArgVar <- resL |> 
        purrr::map(~ attr(.x, "subclass_arg")) |> 
        unlist(use.names = FALSE) |> 
        unique()
    names(resL) <- names(seL)
    res <- dplyr::bind_rows(resL, .id = "Rank") |> 
        dplyr::relocate(.data$Rank, .after = tidyselect::last_col())
    class(res) <- c("lefser_df_clades", class(res))
    attr(res, "pathStrings") <- pathStrings
    attr(res, "tree") <- .toTree(pathStrings)
    attr(res, "lclassf") <- controlVar
    attr(res, "case") <- caseVar
    attr(res, "inputSE") <- seL
    attr(res, "class_arg") <- classArgVar
    attr(res, "subclass_arg") <- subclassArgVar
    return(res)
}

## This functions makes sure that only the taxonomy
## is used for the rowData.
## OTU's or other non-typical taxonomic ranks will not be included.
.getTaxonomyFromPathStr <- function(pathStrings) {
    rgx <- "^k__[^|]+\\|p__[^|]+\\|c__[^|]+\\|o__[^|]+\\|f__[^|]+(\\|g__[^|]+)?(\\|s__[^|]+)?(\\|t__[^|]+)?"
    stringr::str_extract(pathStrings, pattern = rgx)
}

## This function selects pathStrings containing only
## taxa that is differentiallty abundant
.selectPathStrings <- function(se, res) {
    pathStrings <- rownames(se)
    index <- res$features |>
        purrr::map(~ which(stringr::str_detect(pathStrings, .x))) |>
        unlist() |>
        unique() |>
        sort()
    pathStrings <- pathStrings[index]
    return(pathStrings)
}

.rowData2PathStrings <- function(se) {
    xDat <- as.data.frame(SummarizedExperiment::rowData(se))
    pathStrings <- tidyr::unite(
        data = xDat, col = "taxonomy", sep = "|",
        1:tidyselect::last_col()
    ) |> 
        dplyr::pull(.data$taxonomy)
    pathStrings <- sub("(\\|NA)+$", "", pathStrings)
    return(pathStrings)
}

.selectTaxRanks <- function(x) {
    taxNames <- c(
        "kingdom", "phylum", "class", "order", "family",
        "genus", "species", "strain"
    )
    se <- x
    rowDat <- SummarizedExperiment::rowData(se)
    colnames(rowDat) <- stringr::str_to_lower(colnames(rowDat))
    if ("superkingdom" %in% colnames(rowDat)) {
        message("superkingdom --> kingdom")
        colnames(rowDat)[which(colnames(rowDat) == "superkingdom")] <- "kingdom"
    }
    selectCols <- intersect(taxNames, colnames(rowDat))
    SummarizedExperiment::rowData(se) <- rowDat[, selectCols]
    return(se)
}

.appendRankLetter <- function(x) {
    se <- x
    rowDat <- SummarizedExperiment::rowData(se)
    for (i in seq_along(rowDat)) {
        chr_vct <- rowDat[[i]]
        rankName <- colnames(rowDat)[i]
        rankLetter <- dplyr::case_when(
            rankName == "kingdom" ~ "k__",
            rankName == "phylum" ~ "p__",
            rankName == "class" ~ "c__",
            rankName == "order" ~ "o__",
            rankName == "family" ~ "f__",
            rankName == "genus" ~ "g__",
            rankName == "species" ~ "s__",
            rankName == "strain" ~ "t__",
        )
        # rankLetter <- stringr::str_c(
        #     stringr::str_extract(colnames(rowDat)[i], "^\\w"), "__"
        # )
        rowDat[[i]] <- stringr::str_c(rankLetter, chr_vct)
    }
    SummarizedExperiment::rowData(se) <- rowDat
    return(se)
}

.dropFeatures <- function(x) {
    se <- x
    row_data <- as.data.frame(SummarizedExperiment::rowData(se))
    nrow1 <- nrow(row_data)
    for (i in seq_along(row_data)) {
        sumNA <- sum(is.na(row_data[[i]]))
        if (sumNA > 0) {
            message(
                sumNA, " features don't have ", colnames(row_data)[i],
                " information."
            )
        }
    }
    row_data <- tidyr::drop_na(row_data)
    nrow2 <- nrow(row_data)
    if (nrow1 > nrow2) {
        message(
            "Dropped ", nrow1 - nrow2,
            " features without full taxonomy information."
        )
    }
    se <- se[rownames(row_data),]
    pathStrings <- sort(unique(.rowData2PathStrings(se)))
    l1 <- length(pathStrings)
    tips <- stringr::str_extract(pathStrings, "(|\\w__\\w+)?$")
    dupTips <- tips[which(duplicated(tips))]
    pathStrings <- pathStrings[!tips %in% dupTips]
    l2 <- length(pathStrings)
    if (l1 > l2) {
        message(
            "Dropped ", l1 - l2,
            " features with duplicated tip names."
        )
    }
    list(se = se, pathStrings = pathStrings)
}

## Convert a character vector with pathStrings into a cladogram
## These could come from the rownames of a SummarizedExperiment with
## terminal nodes
.toTree <- function(pathStrs) {
    edgeDF <- pathStrs |>
        purrr::map(.pathString2EdgeList) |>
        dplyr::bind_rows() |>
        dplyr::distinct()
    tipLabels <- stringr::str_extract(pathStrs, "[^|]+$")
    nodeLabels <- unique(edgeDF$from)
    idMap <- 1:(length(tipLabels) + length(nodeLabels))
    names(idMap) <- c(tipLabels, nodeLabels)
    edgeMat <- matrix(
        data = c(idMap[edgeDF$from], idMap[edgeDF$to]),
        ncol = 2
    )
    tr <- list(
        edge = edgeMat,
        tip.label = tipLabels,
        node.label = nodeLabels,
        Nnode = length(nodeLabels),
        Ntip = length(tipLabels)
    )
    class(tr) <- "phylo"
    tr
}

## Helper function for .toTree
## Input is a single path string, e.g., "k__bacteria|p_Fusobacteria..."
.pathString2EdgeList <- function(pathStr) {
    pathStrRoot <- stringr::str_c("ROOT|", pathStr)
    chr_vct <- stringr::str_split(pathStrRoot, "\\|")[[1]]
    data.frame(
        from = chr_vct[1:length(chr_vct)-1],
        to = chr_vct[2:length(chr_vct)]
    )
}