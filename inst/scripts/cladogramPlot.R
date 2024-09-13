# library(SummarizedExperiment)

data("zeller14")
z14 <- zeller14[, zeller14$study_condition != "adenoma"]
tn <- get_terminal_nodes(rownames(z14))
z14tn <- z14[tn, ]
# z14tn_ra <- relativeAb(z14tn)

taxDat <- data.frame(x  = rownames(z14tn)) |>
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
rowData(z14tn) <- taxDat

l <- mia::splitByRanks(z14tn)
l <- l[names(l) != "kingdom"]
l <- purrr::map(l, relativeAb)

resl <- purrr::map(l, ~ lefser(relab = .x, groupCol = "study_condition"))
res <- dplyr::bind_rows(resl)

z14tn_ra <- relativeAb(z14tn)
resOriginal <- lefser(relab = z14tn_ra, groupCol = "study_condition")
resOriginal$features  <- stringr::str_extract(resOriginal$features, "[^|]+$")


res2 <- res |> 
    dplyr::filter(!features %in% resOriginal$features) |> 
    dplyr::bind_rows(resOriginal)



pathString <- rownames(z14tn)

index <- res2$features |> 
    purrr::map(~ which(stringr::str_detect(pathString, .x))) |> 
    unlist() |> 
    unique() |> 
    sort()

path_string <- pathString[index]

# res$features <- stringr::str_extract(res$features, "[^|]+$")
res3 <- res2 |> 
    mutate(sample = ifelse(scores > 0, "CRC", "control")) |> 
    mutate(abs = abs(scores)) 

tree <- .toTree(path_string)

labels <- c(tree$tip.label, tree$node.label)

res3$node <- match(res3$features, labels)
res3 <- dplyr::relocate(res3, node)

# library(treeio)

internalNodes <- Ntip(tree) + 1:Nnode(tree)


collapseThem <- map_int(internalNodes, ~ {
    chnods <- offspring(.data = tree, .node = .x, type = "tips")
    if (any(chnods %in% res3$node)) {
        return(NA)
    } else {
        return(.x)
    }
}) |> 
    discard(is.na)

res3 <- res3 |> 
    dplyr::mutate(
        showNodeLabs = dplyr::case_when(
            grepl("[pc]__", features) ~ features, 
            TRUE ~ NA
            # TRUE ~ as.character(node)
        )
    )



# res2 <- dplyr::relocate(res2, node)

x <- rtree(10)

gt <- ggtree(
    tree, layout = "circular",  branch.length = "none", size = 0.2
) %<+% res3 
gt2 <- gt +
    geom_tiplab(aes(label = features), size = 2, geom = "label") +
    geom_tippoint(aes(fill = sample, size = abs), shape = 21) +
    geom_nodepoint(aes(fill = sample, size = abs), shape = 21) +
    geom_text_repel(aes(label = showNodeLabs)) +
    theme(legend.position = "right") + 
    scale_fill_manual(
        values = c("red", "forestgreen"), breaks = c("control", "CRC")
    )

for (i in collapseThem) {
    gt2 <- ggtree::collapse(gt2, node = i)
}
