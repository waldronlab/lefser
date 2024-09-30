
#' LEfSer plot cladogram
#' 
#' \code{lefserPlotClad} plots a cladogram from the results of
#' \code{lefserClades}.
#'
#' @param df An object of class "lefesr_df_clades".
#' @param colors Colors corresponding to class 0 and 1.
#' Options: "c" (colorblind), "l" (lefse), "g" (greyscale).
#' Defaults to "c".
#' This argument also accepts a character(2) with two color names.
#' @param showTipLabels Logical. If TRUE, show tip labels. Default is FALSE.
#' @param showNodeLabels  Label's to be shown in the tree.
#' Options: "p" = phylum, "c" = class, "o" = order,
#' "f" = family, "g" = genus. It can accept several options, e.g., 
#' c("p", "c").
#'
#' @importFrom ggtree %<+%
#'
#' @return A ggtree object.
#' @export
#'
#' @examples
#' data("zeller14")
#' z14 <- zeller14[, zeller14$study_condition != "adenoma"]
#' tn <- get_terminal_nodes(rownames(z14))
#' z14tn <- z14[tn, ]
#' z14tn_ra <- relativeAb(z14tn)
#' resAll <- lefserAllRanks(relab = z14tn_ra, classCol = "study_condition")
#' ggt <- lefserPlotClad(df = resAll)
lefserPlotClad <- function(
        df, colors = "c", showTipLabels = FALSE, showNodeLabels = "p"
) {
    inputClass <- class(df)[1]
    if (inputClass == "lefser_df") {
        message("Working with lefser_df. Consider using lefserAll.")
    } else if (inputClass == "lefser_df_clades") {
        message("Working with lefser_df_clades")
    } else {
        stop(
            "You need an object of class 'lefser_df_class'",
            call. = FALSE
        )
    }
    
    df$features <- .extracTips(df$features)
    
    colors <- .selectPalette(colors)
    tree <- attr(df, "tree")
    controlVar <- attr(df, "lclassf")
    caseVar <- attr(df, "case")
    
    res <- df |>
        dplyr::mutate(
            sample = dplyr::case_when(
                ## This assumes positive values always mean enriched in
                ## the case condition.
                .data[["scores"]] > 0 ~ .env[["caseVar"]],
                TRUE ~ .env[["controlVar"]]
            )
        ) |>
        dplyr::mutate(abs = abs(.data[["scores"]])) |>
        as.data.frame()
    
    labels <- c(tree$tip.label, tree$node.label)
    res$node <- match(res$features, labels)
    dat <- dplyr::relocate(res, node)
    
    internalNodes <- ape::Ntip(tree) + 1:ape::Nnode(tree)
    
    # collapseThem <- purrr::map_int(internalNodes, ~ {
    #     chNods <- treeio::offspring(.data = tree, .node = .x, type = "tips")
    #     if (any(chNods %in% dat$node)) {
    #         return(NA)
    #     } else {
    #         return(.x)
    #     }
    # }) |>
    #     purrr::discard(is.na)
    
    nodLab <- match.arg(
        arg = showNodeLabels,
        choices = c("p", "c", "o", "f", "g"),
        several.ok = TRUE
    )
    nodLabRgx <- paste0("[", paste0(nodLab, collapse = ""), "]__")
    treeData <- dat |>
        dplyr::mutate(
            showNodeLabs = dplyr::case_when(
                grepl(nodLabRgx, features) ~ features,
                TRUE ~ NA
            )
        )
    
    gt <- ggtree::ggtree(
        tree, layout = "circular",  branch.length = "none", size = 0.2
    ) %<+% treeData
    
    gt <- gt +
        ggtree::geom_tippoint(
            mapping = ggtree::aes(fill = sample, size = abs), shape = 21,
            na.rm=TRUE
        ) +
        ggtree::geom_nodepoint(
            mapping = ggtree::aes(fill = sample, size = abs), shape = 21,
            na.rm = TRUE
        )
    
    if (showTipLabels) {
        gt <- gt + 
            ggtree::geom_tiplab(
                mapping = ggtree::aes(label = features), size = 2,
                geom = "text", na.rm=TRUE
            )
    }
    
    gt2 <- gt +
        ggrepel::geom_label_repel(
            mapping = ggtree::aes(label = showNodeLabs),
            na.rm = TRUE
        ) +
        ggtree::scale_fill_manual(
            values = colors, breaks = c(controlVar, caseVar),
            name = "Sample", na.value = NA
        ) +
        ggplot2::scale_size(name = "Absolute\nscore") +
        ggtree::theme(legend.position = "right")
    
    # for (i in collapseThem) {
    #     gt2 <- withCallingHandlers(
    #         ggtree::collapse(gt2, node = i),
    #         warning = function(w) {
    #             if (grepl("collapse", w$message)) {
    #                 invokeRestart("muffleWarning")
    #             }
    #         }
    #     )
    # }
    return(gt2)
}

# Run lefser at all taxonomic levels --------------------------------------

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
#' z14_input <- rowNames2RowData(z14_tn_ra)
#'
#' resAll <- lefserAllRanks(relab = z14tn_input, groupCol = "study_condition")
#' head(resAll)
#' 
lefserClades <- function(relab, ...) {
    se <- .selectTaxRanks(relab)
    se <- .appendRankLetter(se)
    l <- .dropFeatures(se)
    se <- l[["se"]]
    pathStrings <- l[["pathStrings"]]
    seL <- as.list(mia::splitByRanks(se))
    ## Restrict to species. Kingdom would not be informative
    seL <- seL[!names(seL) %in% c("kingdom", "strain")]
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
            "\n>>>> Running lefse at the ", idx, " level.", 
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
    names(resL) <- names(seL)
    res <- dplyr::bind_rows(resL, .id = "Rank") |> 
        dplyr::relocate(.data$Rank, .after = tidyselect::last_col())
    class(res) <- c("lefser_df_clades", class(res))
    attr(res, "pathStrings") <- pathStrings
    attr(res, "tree") <- .toTree(pathStrings)
    attr(res, "lclassf") <- controlVar
    attr(res, "case") <- caseVar
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
    colnames(rowDat)[which(colnames(rowDat) == "superkingdom")] <- "kingdom"
    selectCols <- intersect(taxNames, colnames(rowDat))
    SummarizedExperiment::rowData(se) <- rowDat[, selectCols]
    return(se)
}

.appendRankLetter <- function(x) {
    se <- x
    rowDat <- SummarizedExperiment::rowData(se)
    for (i in seq_along(rowDat)) {
        chr_vct <- rowDat[[i]]
        rankLetter <- stringr::str_c(stringr::str_extract(colnames(rowDat)[i], "^\\w"), "__")
        rowDat[[i]] <- stringr::str_c(rankLetter, chr_vct)
    }
    SummarizedExperiment::rowData(se) <- rowDat
    return(se)
}

.dropFeatures <- function(x) {
    se <- x
    row_data <- as.data.frame(SummarizedExperiment::rowData(se))
    nrow1 <- nrow(row_data)
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

## This function extracts only the last element of the taxonomy
.extracTips <- function(pathStrs) {
    stringr::str_extract(pathStrs, "[^|]+$")
}

#' RowNames to RowData
#' 
#' \code{rowNames2RowData} transforms the taxonomy stored in the row names to 
#' the rowData in a SummarizedExperiment.
#' 
#' @param x A SummarizedExperiment with the features taxonomy in the rownames.
#'
#' @return The same SummarizedExpriment with the taxonomy now in the rowData.
#' @export
#'
#' @examples
#'
#' data("zeller14")
#' 
#' ## Keep only "CRC" and "control" (dichotomous variable)
#' z14 <- zeller14[, zeller14$study_condition %in% c("control", "CRC")]
#' 
#' ## Get terminal nodes
#' tn <- get_terminal_nodes(rownames(z14))
#' z14_tn <- z14[tn, ]
#' 
#' ## Normalize to relative abundance (also known as Total Sum Scaling)
#' z14_tn_ra <- relativeAb(z14_tn)
#' 
#' ## Add the taxonomy to the rowData
#' input_se <- rowNames2RowData(z14_tn_ra)
#' 
rowNames2RowData <- function(x) {
    se <- x
    taxonomy <- .getTaxonomyFromPathStr(rownames(se))
    dataFrame <- data.frame(tax = taxonomy) |>
        tidyr::separate(
            col = "tax", into = paste0("col", 1:10), # Number of taxa is usually seven, so 10 should be more than enough.
            sep = "\\|", extra = "merge", fill = "right"
        ) |>
        purrr::discard(~ all(is.na(.x)))
    ## purrr::map_chr ensures that the a single letter is used per column.
    ## Having two or more letters would trigger and error message from map_chr.
    firstLetter <- purrr::map_chr(dataFrame, ~ {
        taxLvl <- stringr::str_extract(.x, "\\w__")
        unique(taxLvl[which(!is.na(taxLvl))])
    })
    rankNames <- dplyr::case_when(
        firstLetter == "k__" ~ "kingdom",
        firstLetter == "p__" ~ "phylum",
        firstLetter == "c__" ~ "class",
        firstLetter == "o__" ~ "order",
        firstLetter == "f__" ~ "family",
        firstLetter == "g__" ~ "genus",
        firstLetter == "s__" ~ "species",
        firstLetter == "t__" ~ "strain",
    )
    colnames(dataFrame) <- rankNames
    dataFrame <- purrr::modify(dataFrame, ~ sub("^\\w__", "", .x))
    DF <- S4Vectors::DataFrame(dataFrame)
    SummarizedExperiment::rowData(se) <- DF
    return(se)
}
