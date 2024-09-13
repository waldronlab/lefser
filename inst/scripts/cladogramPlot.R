library(SummarizedExperiment)
## Taxonomy is encoded in the rownames
data("zeller14")
z14 <- zeller14[, zeller14$study_condition != "adenoma"]
tn <- get_terminal_nodes(rownames(z14))
z14tn <- z14[tn, ]
z14tn_ra <- relativeAb(z14tn)

res <- lefser(relab = z14tn_ra, groupCol = "study_condition")

taxDat <- data.frame(x  = rownames(z14tn_ra)) |>
    tidyr::separate(
        col = "x", into = paste0("col", 1:10),
        sep = "\\|", extra = "merge", fill = "right"
    ) |>
    purrr::discard(~ all(is.na(.x)))
firstLetter <- purrr::map_chr(taxDat, ~ {
    taxLev <- stringr::str_extract(.x, "\\w__")
    unique(taxLev[which(!is.na(taxLev))])
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

colnames(taxDat) <- rankNames
taxDat$pathString <- rownames(z14tn_ra)

res$pathString <- res$features
res$features <- stringr::str_extract(res$features, "[^|]+$")
res2 <- res |> 
    mutate(sample = ifelse(scores > 0, "CRC", "control")) |> 
    mutate(abs = abs(scores)) 
# eolnames(res2)[which(colnames(res2) == "features")] <- "node"
tree <- .toTree(res2$pathString)
res2$node <- match(res2$features, tree$tip.label)
res2 <- dplyr::relocate(res2, node)

gt <- ggtree(
    tree, layout = "circular",  branch.length = "none", size = 0.25
) %<+% res2 
gt +
    geom_tiplab(aes(label = features), size = 2, geom = "label") +
    geom_tippoint(aes(color = sample, size = abs)) +
    theme(legend.position = "right") + 
    scale_color_manual(
        values = c("purple", "pink"), breaks = c("control", "CRC")
    ) +
    scale_size_continuous(range = c(3, 10))