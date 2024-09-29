suppressMessages(library(lefser))
data(zeller14)
zeller14 <- zeller14[, zeller14$study_condition != "adenoma"]
tn <- get_terminal_nodes(rownames(zeller14))
zeller14tn <- zeller14[tn,]
zeller14tn_ra <- relativeAb(zeller14tn)

res <- lefser(
    zeller14tn_ra, classCol = "study_condition", subclassCol = "age_category"
)
lefserPlot(res)
lefserPlot(res, colors = "l")
lefserPlot(res, colors = "g")
lefserPlot(res, colors = c("purple", "pink"))
sessioninfo::session_info()





# (1) Using classes only
res_class <- lefser(zeller14tn_ra,
                    classCol = "study_condition")
# (2) Using classes and sub-classes
res_subclass <- lefser(zeller14tn_ra,
                    classCol = "study_condition",
                    subclassCol = "age_category")

plot_class <- lefsePlotFeat(res_class, res_class$features[[1]], colors = "g")
plot_subclass <- lefsePlotFeat(res_subclass, res_subclass$features[[2]], colors = "g")

plot_class <- lefsePlotFeat(res_class, res_class$features[[1]], colors = c("orange", "blue"))
plot_subclass <- lefsePlotFeat(res_subclass, res_subclass$features[[1]], colors = c("orange", "blue"))
plot_subclass

lefserPlot(res_subclass, colors = "")






