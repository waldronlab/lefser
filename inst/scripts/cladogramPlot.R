# suppressPackageStartupMessages(library(lefser))
data("zeller14")
z14 <- zeller14[, zeller14$study_condition != "adenoma"]
tn <- get_terminal_nodes(rownames(z14))
z14tn <- z14[tn, ]
z14tn_ra <- relativeAb(z14tn)

res <- lefser(z14tn_ra, groupCol = "study_condition")
resAll <- lefserAllRanks(relab = z14tn_ra, groupCol = "study_condition")

x <- lefserPlotClad(df = res)
x
y <- lefserPlotClad(df = resAll)
y
sessioninfo::session_info()
