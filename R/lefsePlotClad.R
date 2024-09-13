lefsePlotClad <- function(x) {
    x
}

.pathString2DataFrame <- function(pathString) {
    s <- stringr::str_c("ROOT|", pathString)
    chr_vct <- stringr::str_split(s, "\\|")[[1]]
    data.frame(
        from = chr_vct[1:length(chr_vct)-1],
        to = chr_vct[2:length(chr_vct)]
    )
}

.toTree <- function(v) {
    edgeDF <- v |> 
        purrr::map(.pathString2DataFrame) |> 
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
    # treeio::as.treedata(tr)
}