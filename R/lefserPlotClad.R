
# Functions for plotting a cladogram --------------------------------------

#' LEfSer plot cladogram
#' 
#' \code{lefserPlotClad} plots a cladogram from the results of 
#' `lefser` or `lefserAllRanks`
#'
#' @param df An object of class "lefser_df" or "lefesr_df_all".
#' @param colors Colors corresponding to class 0 and 1.
#' Options: "c" (colorblind), "l" (lefse), "g" (greyscale).
#' Defaults to "c". This argument also accepts a character(2) with two color names.
#' @param showTipLabels Logical. If TRUE, show tip labels. Default is FALSE.
#' @param showNodeLabels Options: "p" = phylum, "c" = class, "o" = order,
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
#' resAll <- lefserAllRanks(relab = z14tn_ra, groupCol = "study_condition")
#' ggt <- lefserPlotClad(df = resAll)
lefserPlotClad <- function(
        df, colors = "c", showTipLabels = FALSE, showNodeLabels = "p"
) {
    inputClass <- class(df)[1]
    if (inputClass == "lefser_df") {
        message("Woriking with lefser_df. Consider using lefserAll.")
        # df$features <- .extracTips(df$features)
    } else if (inputClass == "lefser_df_all") {
        message("Working with lefser_df_all")
        ## .extractTips should be use here as well
        ## The feature names format should use full taxonomy
    } else {
        stop(
            "You need an object of class 'lefser_df_all'",
            call. = FALSE
        )
    }
    
    df$features <- .extracTips(df$features)
    
    colors <- .selectPalette(colors)
    tree <- attr(df, "tree")
    controlVar <- attr(df, "lgroupf")
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
    collapseThem <- purrr::map_int(internalNodes, ~ {
        chNods <- treeio::offspring(.data = tree, .node = .x, type = "tips")
        if (any(chNods %in% dat$node)) {
            return(NA)
        } else {
            return(.x)
        }
    }) |>
        purrr::discard(is.na)
    
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
    # return(treeData)
    
    gt <- ggtree::ggtree(
        tree, layout = "circular",  branch.length = "none", size = 0.2
    ) %<+% treeData
    
    if (showTipLabels) {
        gt <- gt + 
            ggtree::geom_tiplab(
                mapping = ggtree::aes(label = features), size = 2,
                geom = "text", na.rm=TRUE
            )
    }
    
    gt2 <- gt +
        ggtree::geom_tippoint(
            mapping = ggtree::aes(fill = sample, size = abs), shape = 21,
            na.rm=TRUE
        ) +
        ggtree::geom_nodepoint(
            mapping = ggtree::aes(fill = sample, size = abs), shape = 21,
            na.rm = TRUE
        ) +
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
    for (i in collapseThem) {
        gt2 <- withCallingHandlers(
            ggtree::collapse(gt2, node = i),
            warning = function(w) {
                if (grepl("collapse", w$message)) {
                    invokeRestart("muffleWarning")
                }
            }
        )
    }
    return(gt2)
}

# Run lefser at all taxonomic levels --------------------------------------

#' Run lefser on all taxonomic levels
#'
#' @param relab A SummarizedExperiment.
#' @param ... Arguments passed to the \code{lefser} function.
#'
#' @return An object of class 'lefser_df_all' and 'data.frame'.
#' @export
#'
#' @examples
#' 
#' data("zeller14")
#' z14 <- zeller14[, zeller14$study_condition != "adenoma"]
#' tn <- get_terminal_nodes(rownames(z14))
#' z14tn <- z14[tn, ]
#' z14tn_ra <- relativeAb(z14tn)
#'
#' resAll <- lefserAllRanks(relab = z14tn_ra, groupCol = "study_condition")
#'
lefserAllRanks <- function(relab,...) {
    ## Feature names should have the full taxonomy
    se <- .rowNames2RowData(relab)
    seL <- mia::splitByRanks(se)
    ## The kingdom level is not needed
    ## The mia package doesn't support strain.
    seL <- seL[names(seL) != "kingdom"]
    seL <- purrr::map(seL, ~ {
        seVar <- .x
        rowDat <- as.data.frame(SummarizedExperiment::rowData(seVar))
        rowDat <- purrr::discard(rowDat, function(x) all(is.na(x)))
        rowDat <- S4Vectors::DataFrame(rowDat)
        SummarizedExperiment::rowData(seVar) <- rowDat
        seVar
    })
    for (i in seq_along(seL)) {
        rownames(seL[[i]]) <- .lognRowNames(seL[[i]])
    }
    
    res <- seL |>
        purrr::map(function(x, ...) lefser(relab = x,...), ...) |>
        dplyr::bind_rows()
    resOriginal <- lefser(relab, ...)
    ## Get only tip names (full names with full taxonomy are too long).
    # resOriginal$features  <- stringr::str_extract(
    #     resOriginal$features, "[^|]+$"
    # )
    res <- res |>
        ## Avoid repeating features.
        dplyr::filter(!.data[["features"]] %in% resOriginal$features) |>
        ## Features not supported by mia are added (strain, OTUs, etc.)
        dplyr::bind_rows(resOriginal)

    controlVar <- attr(resOriginal, "lgroupf")
    caseVar <- attr(resOriginal, "case")

    class(res) <- c("lefser_df_all", class(res))

    ## These pathStrings could be used in the plotting function instead (or not)
    pathStrings <- .selectPathStrings(relab, res)
    attr(res, "pathStrings") <- pathStrings
    attr(res, "tree") <- .toTree(pathStrings)

    attr(res, "lgroupf") <- controlVar
    attr(res, "case") <- caseVar
    return(res)
}

## Add taxonomic information to rowData
## This step is necessary for mia to work
.rowNames2RowData <- function(x) {
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
    DF <- S4Vectors::DataFrame(dataFrame)
    SummarizedExperiment::rowData(se) <- DF
    return(se)
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

# Create cladogram --------------------------------------------------------
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


# Utils -------------------------------------------------------------------
.lognRowNames <- function(se) {
    dat <- SummarizedExperiment::rowData(se) |> 
        as.data.frame() |> 
        tibble::rownames_to_column(var = "rowname") |> 
        dplyr::relocate(.data[["rowname"]])
    lastColLgl <- all(dat[[colnames(dat)[ncol(dat)]]] == dat[["rowname"]])
    if (lastColLgl) {
        dat <- dat[, which(colnames(dat) != "rowname")]
        output <- dat |> 
            tidyr::unite(
                col = "features", 1:tidyselect::last_col(), 
                sep = "|", remove = TRUE, 
            ) |> 
            dplyr::pull(.data[["features"]])
    } 
    return(output)
}
