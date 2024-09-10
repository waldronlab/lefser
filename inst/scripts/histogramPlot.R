suppressMessages(library(lefser))
data(zeller14)
zeller14 <- zeller14[, zeller14$study_condition != "adenoma"]
tn <- get_terminal_nodes(rownames(zeller14))
zeller14tn <- zeller14[tn,]
zeller14tn_ra <- relativeAb(zeller14tn)

# (1) Using classes only
res_group <- lefser(zeller14tn_ra,
                    groupCol = "study_condition")
# (2) Using classes and sub-classes
res_block <- lefser(zeller14tn_ra,
                    groupCol = "study_condition",
                    blockCol = "age_category")
plot_group <- lefsePlotFeat(res_group, res_group$features[[1]])
plot_group
plot_block <- lefsePlotFeat(res_block, res_block$features[[2]])
plot_block
