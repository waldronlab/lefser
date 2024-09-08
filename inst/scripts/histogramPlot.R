library(lefser)
library(ggplot2)
library(dplyr)
data(zeller14)
zeller14 <- zeller14[, zeller14$study_condition != "adenoma"]
tn <- get_terminal_nodes(rownames(zeller14))
zeller14tn <- zeller14[tn,]
zeller14tn_ra <- relativeAb(zeller14tn)

res_block <- lefser(zeller14tn_ra, 
                    groupCol = "study_condition", 
                    blockCol = "age_category")

dat <- lefserPlotHistogram(
    res_block, zeller14tn_ra, groupCol = "study_condition",
    blockCol = "age_category"
)

bacName <- "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Ruminococcus|s__Ruminococcus_sp_5_1_39BFAA|t__GCF_000159975"
dat |> 
    filter(features == bacName) |> 
    arrange(desc(study_condition), desc(age_category), desc(abundance)) |> 
    mutate(sample = forcats::fct_inorder(sample)) |> 
    ggplot(aes(sample, abundance)) +
    geom_col(aes(fill = age_category)) +
    facet_wrap(~ study_condition, scales = "free_x") +
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
