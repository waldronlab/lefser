
## Taxonomy is encoded in the rownames
data("zeller14")
z14 <- zeller14
tn <- get_terminal_nodes(rownames(z14))
z14tn <- z14[tn,]
z14tn_ra <- relativeAb(z14tn)

## Let's transfer this taxonomy to the rowData
dat <- data.frame(x  = rownames(z14tn_ra)) |>
    tidyr::separate(
        col = "x", into = paste0("col", 1:10),
        sep = "\\|", extra = "merge", fill = "right"
    ) |>
    purrr::discard(~ all(is.na(.x)))
firstLetter <- purrr::map_chr(dat, ~ {
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
        firstLetter == "t__" ~ "strain"
    )

colnames(dat) <- rankNames
dat$pathString <- rownames(z14tn_ra)

#
# tr <- toTree(dat)
# data.tree::FromDataFrameTable(dat)
# data_tree <- as.Node(dat, pathName = "pathString", pathDelimiter = "|")
# nw <- data.tree::ToNewick(data_tree)
# sl <- data.tree::ToListSimple(data_tree)

t1 <- taxa2phylo(paste0("ROOT|", dat$pathString))

ggtree(t1, layout = "circular")













taxonname2edgelist <- function(taxon) {
    # print(taxon)
    v = strsplit(taxon,'\\|')[[1]]
    v = v[!v=='NA']
    if(length(v)>1) {
        lv = length(v)
        df = data.frame(from=v[1:(lv-1)], to = v[2:lv])
    } else {
        df = data.frame()
    }
    df
}


taxa2edgelist <- function(taxa) {
    taxa_edgelist <- lapply(taxa,taxonname2edgelist)
    df = unique(do.call(rbind, taxa_edgelist))
    return(df)
}

taxa2phylo <- function(taxa) {
    edgelist = taxa2edgelist(taxa)
    edgelist = as.matrix(edgelist)

    edgelist = edgelist[!is.na(edgelist[,1]) & !is.na(edgelist[,2]),]

    from <- edgelist[,1]
    to <- edgelist[,2]
    ids <- unique(c(edgelist[,1], edgelist[,2]))

    tip.label <- setdiff(ids, from)
    node.label <- unique(from)

    # make a map from taxonomy ID to internal 1:n ids
    idmap <- 1:(length(tip.label) + length(node.label))
    names(idmap) <- c(tip.label, node.label)

    # make a phylo object
    tree <- list(
        edge       = matrix(c(idmap[as.character(from)], idmap[as.character(to)]), ncol=2),
        tip.label  = unname(tip.label),
        node.label = unname(node.label),
        Nnode      = length(node.label)
    )
    class(tree) <- 'phylo'

    tree

}

t1 <- taxa2phylo(nms)



