utils::globalVariables(c("Names", "scores"))

#' Plots results from `lefser` function
#'
#' `lefserPlot` function displays effect sizes for differentially expressed microorganisms
#' and whether they are more abundant in '0' or '1' sample group.
#'
#' @param df
#' Data frame produced by `lefser`.
#'
#' @param colors character(2) The two colors corresponding to class 0 and 1,
#' respectively. Defaults to `c("red", "forestgreen")`.
#' 
#' @param trim.names If `TRUE` extracts the most specific taxonomic rank of organism.
#' 
#' @return
#' Function returns plot of effect size scores produed by `lefser`.
#' Positive scores represent microorganisms with that are more abundant in class '1'.
#' Negative scores represent microorganisms with that are more abundant in class '0'.
#'
#' @export
#' @importFrom ggplot2 ggplot aes ylab theme element_blank element_text
#' @importFrom ggplot2 geom_bar coord_flip scale_fill_manual
#'
#' @examples
#' example("lefser")
#' lefserPlot(res_group)
lefserPlot <- function(df, colors = c("red", "forestgreen"), 
                       trim.names = TRUE) {
  df <- .trunc(df, trim.names)
  group <- ifelse(df$scores > 0, "treatment", "control")
  df$group <- as.factor(group)
  plt <-
    ggplot(df, aes(reorder(Names, scores), scores)) + ylab("LDA SCORE (log 10)") +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 11, face = "bold"),
      axis.text.y  = element_text(
        vjust = 0.7,
        size = 9,
        face = "bold"
      ),
      axis.text.x  = element_text(
        vjust = 0.7,
        size = 9,
        face = "bold"
      ),
      plot.title = element_text(
        hjust = 0.5,
        size = 13,
        face = "bold"
      )
    ) +
    geom_bar(stat = "identity", aes(fill = group)) +
    scale_fill_manual(values = colors) +
    coord_flip()
  return(plt)
}
