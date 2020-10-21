library(lefser)

zell <- zeller14[, zeller14$study_condition != "adenoma"]
exprs <- assay(zell, i = 1L)
datas <- data.frame(
    condition = lefser:::.numeric01(zell$study_condition),
    age = lefser:::.numeric01(zell$age_category),
    ID = colnames(exprs),
    t(exprs)
)

text <- mapply(function(x, y) {
    paste0(x, "\t", y)
    }, x = names(datas), y = sapply(datas, paste0, collapse = "\t")
)

## file for http://huttenhower.sph.harvard.edu/galaxy/ LEfSe
writeLines(text, con = file("~/data/lefser/setupdata.txt"))

## compared to results from lefser
res <- lefser(zell, groupCol = "study_condition", blockCol = "age_category")
lefserPlot(res)
