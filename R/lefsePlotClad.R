
# Functions for plotting a cladogram --------------------------------------

#' Plot Cladogram
#'
#' @param x An object of class "lefser_df" or "lefesr_df_all".
#'
#' @importFrom ggtree %<+%
#'
#' @return A ggtree object.
#' @export
#'
lefsePlotClad <- function(x) {
    # if (!"lefser_df_all" %in% class(x)) {
    #     stop(
    #         "You need an object of class 'lefser_df_all'",
    #         call. = FALSE
    #     )
    # }
    tree <- attr(x, "tree")
    control <- attr(x, "lgroupf")
    case <- attr(x, "case")

    labels <- c(tree$tip.label, tree$node.label)
    x$node <- match(x$features, labels)
    dat <- dplyr::relocate(x, node)
    internalNodes <- ape::Ntip(tree) + 1:ape::Nnode(tree)
    collapseThem <- purrr::map_int(internalNodes, ~ {
        chNods <- treeio::offspring(.data = tree, .node = .x, type = "tips")
        if (any(chNods %in% x$node)) {
            return(NA)
        } else {
            return(.x)
        }
    }) |>
        purrr::discard(is.na)

    treeData <- dat |>
        dplyr::mutate(
            showNodeLabs = dplyr::case_when(
                grepl("[p]__", features) ~ features,
                TRUE ~ NA
            )
        )
    gt <- ggtree::ggtree(
        tree, layout = "circular",  branch.length = "none", size = 0.2
    ) %<+% treeData
    gt2 <- gt +
        ggtree::geom_tiplab(
            mapping = ggtree::aes(label = features), size = 2, geom = "text"
        ) +
        ggtree::geom_tippoint(
            mapping = ggtree::aes(fill = sample, size = abs), shape = 21
        ) +
        ggtree::geom_nodepoint(
            mapping = ggtree::aes(fill = sample, size = abs), shape = 21
        ) +
        ggrepel::geom_label_repel(
            mapping = ggtree::aes(label = showNodeLabs)
        ) +
        ggtree::theme(legend.position = "right") +
        ggtree::scale_fill_manual(
            values = c("red", "forestgreen"), breaks = c(control, case)
        )

    for (i in collapseThem) {
        gt2 <- ggtree::collapse(gt2, node = i)
    }
    gt2
}


.cladPlotTips <- function() {

}

.cladPlotAll <- function() {

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
lefserAllRanks <- function(relab,...) {
    se <- .rowNames2RowData(relab)
    seL <- mia::splitByRanks(se)
    ## The kingdom level is not needed
    ## The mia package doesn't support strain.
    seL <- seL[names(seL) != "kingdom"]
    res <- seL |>
        purrr::map(function(x, ...) lefser(relab = x,...), ...) |>
        dplyr::bind_rows()
    resOriginal <- lefser(relab, ...)
    ## Get only tip names (full names with full taxonomy are too long).
    resOriginal$features  <- stringr::str_extract(
        resOriginal$features, "[^|]+$"
    )
    res <- res |>
        ## Avoid repeating features.
        dplyr::filter(!.data[["features"]] %in% resOriginal$features) |>
        ## Features not supported by mia are added (strain, OTUs, etc.)
        dplyr::bind_rows(resOriginal)

    controlVar <- attr(resOriginal, "lgroupf")
    caseVar <- attr(resOriginal, "case")

    ## This code could be added to the plotting function
    ## To avoid modfying the current lefser function
    resOutput <- res |>
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

    class(resOutput) <- c("lefser_df_all", class(resOutput))

    ## These pathStrings could be used in the plotting function instead (or not)
    pathStrings <- .selectPathStrings(relab, res)
    attr(resOutput, "pathStrings") <- pathStrings
    attr(resOutput, "tree") <- .toTree(pathStrings)

    attr(resOutput, "lgroupf") <- controlVar
    attr(resOutput, "case") <- caseVar
    return(resOutput)
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
