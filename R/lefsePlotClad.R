#' Run lefser on all taxonomic levels
#'
#' @param relab A SummarizedExperiment.
#' @param ... Arguments passed to the \code{lefser} function.
#'
#' @return An object of class 'lefser_df_all' and 'data.frame'.
#' @export
#'
lefserAllRanks <- function(relab,...) {
    se <- rowNames2RowData(relab)
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
    pathStrings <- .selectPathStrings(relab, res)
    
    resOutput <- res |>
        dplyr::mutate(
            sample = dplyr::case_when(
                ## This assumes positive values always mean enriched in
                ## the case condition.
                .data[["scores"]] > 0 ~ .env[["case"]],
                TRUE ~ .env[["controlVar"]]
            )
        ) |> 
        dplyr::mutate(abs = abs(.data[["scores"]])) |> 
        as.data.frame()
    
    class(resOutput) <- c("lefser_df_all", class(resOutput))
    attr(resOutput, "pathStrings") <- pathStrings
    attr(resOutput, "tree") <- .toTree(pathStrings)
    attr(resOutput, "lgroupf") <- control
    attr(resOutput, "case") <- caseVar
    
    return(resOutput)
}

#' Plot Cladogram
#'
#' @param x An object of class "lefesr_df_all".
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
    browser()
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

## Add taxonomic information to rowData
rowNames2RowData <- function(x) {
    se <- x
    dataFrame <- data.frame("rn" = rownames(se)) |> 
        tidyr::separate(
            col = "rn", into = paste0("col", 1:10),
            sep = "\\|", extra = "merge", fill = "right"
        ) |> 
        purrr::discard(~ all(is.na(.x)))
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
    SummarizedExperiment::rowData(se) <- dataFrame
    return(se)
}

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

## Convert a character vector into a cladogram based on taxonomy
.toTree <- function(v) {
    edgeDF <- v |> 
        purrr::map(.pathString2EdgeList) |> 
        dplyr::bind_rows() |> 
        dplyr::distinct()
    tipLabels <- stringr::str_extract(v, "[^|]+$")
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
.pathString2EdgeList <- function(pathString) {
    s <- stringr::str_c("ROOT|", pathString)
    chr_vct <- stringr::str_split(s, "\\|")[[1]]
    data.frame(
        from = chr_vct[1:length(chr_vct)-1],
        to = chr_vct[2:length(chr_vct)]
    )
}
