#' Utility function to calculate relative abundances
#'
#' @description
#' The function calculates the column totals and divides each value within the
#' column by the respective column total.
#'
#' @inheritParams lefser
#'
#' @param se A SummarizedExperiment object with counts
#'
#' @examples
#'
#' se <- SummarizedExperiment(
#'     assays = list(
#'         counts = matrix(
#'             rep(1, 4), ncol = 1, dimnames = list(LETTERS[1:4], "SAMP")
#'         )
#'     )
#' )
#' assay(se)
#' assay(relativeAb(se))
#'
#' @description
#' This function calculates the relative abundance of each feature in the SummarizedExperiment
#' object containing count data, expressed as counts per million (CPM)
#'
#' @returns returns a new SummarizedExperiment object with counts per million
#' calculated and added as a new assay named rel_abs.
#'
#' @export
relativeAb <- function(se, assay = 1L) {
  assay_data <- assay(se, i = assay)
  csums <- colSums(assay_data)
  div <- matrix(rep(csums, each = nrow(assay_data)), ncol = ncol(assay_data))
  res <- assay_data / div * 1e6
  assaylist <- assays(se)
  newalist <- append(
    assaylist, values = S4Vectors::SimpleList(rel_abs = res), after = 0L
  )
  assays(se) <- newalist
  se
}

#' Identify which elements of a string are terminal nodes
#'
#' @param string A character vector of strings to check for terminal nodes
#'
#' @return A logical vector indicating which elements of the string are terminal
#' nodes
#' @description
#' A terminal node in a taxonomy does not have any child nodes. For example, a
#' species is a terminal node if there are no subspecies or strains that belong
#' to that species. This function identifies which elements of a vector are terminal
#' nodes simply by checking whether that element appears as a substring in any other
#' element of the vector.
#'
#' @export
#'
#' @examples
#' # What does it do?
#' data("zeller14")
#' rownames(zeller14)[988:989]
#' get_terminal_nodes(rownames(zeller14)[988:989])
#' # How do I use it to keep only terminal nodes for a lefser analysis?
#' terminal_nodes <- get_terminal_nodes(rownames(zeller14))
#' zeller14sub <- zeller14[terminal_nodes, ]
#' # Then continue with your analysis!
get_terminal_nodes <- function(string) {
  terminal_nodes <- logical(length(string)) # Initialize logical vector
  for (i in seq_along(string)) {
    # Check if the string appears as a substring in any other strings
    if (!any(grepl(string[i], string[-i], fixed = TRUE))) {
      terminal_nodes[i] <- TRUE
    }
  }
  return(terminal_nodes)
}


# Truncate the feature name
.trunc <- function(scores_df, trim.names){
    Names <- gsub("`", "", scores_df[["features"]])
    if (trim.names) {
        listNames <- strsplit(Names, "\\||\\.")
        Names <- vapply(listNames, tail, character(1L), 1L)
    }
    scores_df[["features"]] <- Names
    return(scores_df)
}

## A function for selecting colors for the plots
.selectPalette <- function(x = "colorblind") {
    if (is.character(x) && length(x) == 1) {
        sel <- match.arg(x, choices = c("colorblind", "lefse", "greyscale"))
        message("Using palette: ", sel)
        presetColors <- list(
            colorblind = c("#E57A77", "#7CA1CC"),
            lefse = c("red", "forestgreen"),
            greyscale = c("grey30", "grey60")
        )
        selectedColor <- presetColors[[sel]]
        return(selectedColor)
    } else {
        return(x)
    }
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



