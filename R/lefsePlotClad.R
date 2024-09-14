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
    seL <- seL[names(seL) != "kingdom"]
    res <- seL |> 
        purrr::map(function(x, ...) lefser(relab = x,...), ...) |> 
        dplyr::bind_rows()
    resOriginal <- lefser(relab, ...)
    resOriginal$features  <- stringr::str_extract(resOriginal$features, "[^|]+$")
    
    res2 <- res |> 
        dplyr::filter(!.data$features %in% resOriginal$features) |> 
        dplyr::bind_rows(resOriginal)
    
    grp <- attr(resOriginal, "groups")
    control <- grp[grp == attr(resOriginal, "lgroupf")]
    case <- grp[grp != attr(resOriginal, "lgroupf")]
    
    pathString <- rownames(relab)
    index <- res2$features |> 
        purrr::map(~ which(stringr::str_detect(pathString, .x))) |> 
        unlist() |> 
        unique() |> 
        sort()
    pathString <- pathString[index]
    
    res3 <- res2 |> 
        dplyr::mutate(sample = ifelse(.data$scores > 0, .env$case, .env$control)) |> 
        dplyr::mutate(abs = abs(scores)) |> 
        as.data.frame()
    class(res3) <- c("lefser_df_all", class(res3))
    
    attr(res3, "pathString") <- pathString
    attr(res3, "tree") <- .toTree(pathString)
    attr(res3, "case") <- case
    attr(res3, "control") <- control
    
    return(res3)
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
    
    if (!"lefser_df_all" %in% class(x)) {
        stop(
            "You need an object of class 'lefser_df_all'",
            call. = FALSE
        )
    }
    
    tree <- attr(x, "tree")
    control <- attr(x, "control")
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
        # ggtree::geom_label(aes(label = showNodeLabs)) +
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

## Convert a character vector into a 
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
