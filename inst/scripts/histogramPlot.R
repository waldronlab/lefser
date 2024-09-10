suppressMessages(library(lefser))
data(zeller14)
zeller14 <- zeller14[, zeller14$study_condition != "adenoma"]
tn <- get_terminal_nodes(rownames(zeller14))
zeller14tn <- zeller14[tn,]
zeller14tn_ra <- relativeAb(zeller14tn)

res <- lefser(
    zeller14tn_ra, groupCol = "study_condition", blockCol = "age_category"
)
lefserPlot(res)
lefserPlot(res, colors = "l")
lefserPlot(res, colors = "g")
lefserPlot(res, colors = c("purple", "pink"))
sessioninfo::session_info()





# (1) Using classes only
res_group <- lefser(zeller14tn_ra,
                    groupCol = "study_condition")
# (2) Using classes and sub-classes
res_block <- lefser(zeller14tn_ra,
                    groupCol = "study_condition",
                    blockCol = "age_category")

plot_group <- lefsePlotFeat(res_group, res_group$features[[1]], colors = "g")
plot_block <- lefsePlotFeat(res_block, res_block$features[[2]], colors = "g")

plot_group <- lefsePlotFeat(res_group, res_group$features[[1]], colors = c("orange", "blue"))
plot_block <- lefsePlotFeat(res_block, res_block$features[[1]], colors = c("orange", "blue"))
plot_block

lefserPlot(res_block, colors = "")






