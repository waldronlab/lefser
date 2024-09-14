suppressMessages(library(lefser))
data(zeller14)
zeller14 <- zeller14[, zeller14$study_condition != "adenoma"]
tn <- get_terminal_nodes(rownames(zeller14))
zeller14tn <- zeller14[tn,]
zeller14tn_ra <- relativeAb(zeller14tn)

# (1) Using classes only
res_class <- lefser(zeller14tn_ra,
                    classCol = "study_condition")
# (2) Using classes and sub-classes
res_subclass <- lefser(zeller14tn_ra,
                    classCol = "study_condition",
                    subclassCol = "age_category")
plot_class <- lefsePlotFeat(res_class, res_class$features[[1]])
plot_class
plot_subclass <- lefsePlotFeat(res_subclass, res_subclass$features[[2]])
plot_subclass
