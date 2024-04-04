utils::globalVariables(c("Names", "scores"))

#' Plots results from `lefser` function
#'
#' This function plots the biomarkers found by LEfSe, that are ranked according
#' to their effect sizes and linked to their abundance in each class.
#'
#' @import ggplot2
#' @importFrom dplyr %>% arrange mutate
#' @importFrom utils head tail
#'
#' @param df Data frame produced by \code{lefser}.
#' @param colors A character(2). Colors corresponding to class 0 and 1. 
#' Defaults to `c("red", "forestgreen")`.
#' @param trim.names Under the default (`TRUE`), this function extracts the 
#' most specific taxonomic rank of organism.
#' @param label.font.size A numeric(1). The font size of the feature labels.
#' The default is `3`.
#'
#' @return
#' Function returns plot of effect size scores produced by \code{lefser}.
#' Positive scores represent the biomarker is more abundant in class '1'.
#' Negative scores represent the biomarker is more abundant in class '0'.
#'
#' @examples
#' example("lefser")
#' lefserPlot(res_group)
#' 
#' @export
lefserPlot <- function(df, 
                       colors = c("red", "forestgreen"),
                       trim.names = TRUE,
                       label.font.size = 3) {
    
    df <- .trunc(df, trim.names)
    groups <- attr(df, "groups")
  
    ## Create the `group` column
     if (!is.null(groups)) {
         group <- ifelse(df$scores > 0, tail(groups, 1), head(groups, 1))
         df$group <- factor(group, levels = groups)
    } else {
        group <- ifelse(df$scores > 0, 1, 0)
        df$group <- as.factor(group)
    }
  
    ## Add the `order` column based on the scores 
    ## To make duplicated features behave independently
    df <- df %>%
        arrange(scores) %>%
        mutate(order = seq_len(nrow(.)))

    plt <-
        ggplot(df, aes(factor(order), scores, width = 0.75)) + # Plot same x-axis values separately
        ylab("LDA SCORE (log 10)") +
        theme_void() +
        theme(
            axis.title.y = element_blank(),
            axis.title.x = element_text(size = 11),
            axis.text.y  = element_blank(),
            axis.text.x  = element_text(vjust = 0.7, size = 9)) +
        geom_bar(
            stat = "identity", aes(fill = group), color = "black", size = 0.3) +
        theme(    # Legends
            legend.position = "top",
            legend.title = element_blank(),
            legend.key.height = unit(0.07, 'cm'),
            legend.key.width = unit(0.6, 'cm')) +
        scale_fill_manual(values = colors) +
        geom_text(    # Feature labeling
            aes(y = 0, label = Names), 
            hjust = ifelse(df$scores < 0, 0, 1),
            nudge_y = ifelse(df$scores < 0, 0.1, -0.1),
            color = "black",
            size = label.font.size) +    
        theme(    # Guide lines
            panel.grid.major.x = element_line(
                color = "grey", size = 0.5, linetype = "dotted"),
            panel.grid.minor.x = element_line(
                color = "grey", size = 0.5, linetype = "dotted")) +
        coord_flip()
      
      return(plt)
}
