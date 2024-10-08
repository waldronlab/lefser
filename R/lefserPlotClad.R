
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
#' "f" = family, "g" = genus, "s" = species, "t" = strain.
#' It can accept several options, e.g., c("p", "c").
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
#' resCl <- lefserClades(relab = z14tn_ra, classCol = "study_condition")
#' ggt <- lefserPlotClad(df = resCl)
lefserPlotClad <- function(
        df, colors = "c", showTipLabels = FALSE, showNodeLabels = "p"
) {
    inputClass <- class(df)[1]
    if (inputClass != "lefser_df_clades") {
        stop(
            "An object of class 'lefser_df_class' is needed.",
            " Please use the `lefserClades` function.",
            call. = FALSE
        )
    }
    
    df$features <- .extractTips(df$features)
    
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
        choices = c("p", "c", "o", "f", "g", "s", "t"),
        several.ok = TRUE
    )
    nodLabRgx <- paste0("[", paste0(nodLab, collapse = ""), "]__")
    treeData <- dat |>
        dplyr::mutate(abs = round(.data$abs, 1)) |> 
        dplyr::mutate(
            showNodeLabs = dplyr::case_when(
                grepl(nodLabRgx, features) ~ features,
                TRUE ~ NA
            )
        )
    minval <- floor(min(treeData$abs, na.rm = TRUE))
    maxval <- floor(max(treeData$abs, na.rm = TRUE))
    
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
        ggplot2::scale_size(
            name = "Absolute\nLDA score",
            limits = c(minval, maxval),
            # range = range(seq_along(seq(minval, maxval, 1))),
            range = c(minval, maxval),
            breaks = seq(minval, maxval, 1)
            
        ) +
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

## This function extracts only the last element of the taxonomy
.extractTips <- function(pathStrs) {
    stringr::str_extract(pathStrs, "[^|]+$")
}
