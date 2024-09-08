# library(lefser)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(purrr)
data(zeller14)
zeller14 <- zeller14[, zeller14$study_condition != "adenoma"]
tn <- get_terminal_nodes(rownames(zeller14))
zeller14tn <- zeller14[tn,]
zeller14tn_ra <- relativeAb(zeller14tn)

res <- lefser(
    zeller14tn_ra, 
    groupCol = "study_condition",
    blockCol = "age_category"
)

dat <- .prepareDataHistogram(res)

dat2 <- dat |> 
    filter(features == bacName) |> 
    arrange(desc(study_condition), desc(age_category), desc(abundance)) |> 
    mutate(sample = forcats::fct_inorder(sample)) |> 
    mutate(classes = paste0(study_condition, "--", age_category))
myClasses <- unique(dat2$classes)

vals <- map(myClasses, ~ {
    subDat2 <- dat2 |> 
        filter(classes == .x)
    sumDat <- summarise(subDat2, mean = mean(abundance))
    meanValue <- pull(sumDat, mean)
    sampleValue <- pull(subDat2, sample)
    min = head(sampleValue, 1)
    max = tail(sampleValue, 1)
    data.frame(x1 = min, x2 = max, y1 = meanValue, y2 = meanValue)
} )



x <- split(dat2, dat2$classes) |> 
    map(~ pull(.x, sample))

unlist(vals[[2]][1:2]) %in% x[[1]]


dat2 |> 
    filter(features == bacName) |> 
    ggplot(aes(sample, abundance)) +
    geom_col(aes(fill = age_category), width = 1) +
    geom_segment(
        data = vals[[1]], mapping = aes(x = x1, xend = x2, y = y1, yend = y2)
    ) +
    geom_segment(
        data = vals[[2]], mapping = aes(x = x1, xend = x2, y = y1, yend = y2)
    ) +
    geom_segment(
        data = vals[[3]], mapping = aes(x = x1, xend = x2, y = y1, yend = y2)
    ) +
    geom_segment(
        data = vals[[4]], mapping = aes(x = x1, xend = x2, y = y1, yend = y2),
    ) +
    annotate(
        "label", x = vals[[1]]$x1, y=max(dat2$abundance), label = "Control",
        hjust = 0, vjust = 1 
    ) +
    annotate(
        "label", x = vals[[3]]$x1, y=max(dat2$abundance), label = "CRC",
        hjust = 0, vjust = 1
    ) +
    # facet_wrap(~ study_condition, scales = "free_x", nrow = 1) +
    theme_bw() +
    scale_y_continuous(
        labels = scales::scientific, expand = c(0, 0)
    ) +
    scale_fill_manual(
        values = c("#189aa7", "#c16973")
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank()
    )

