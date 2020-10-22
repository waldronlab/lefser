library(lefser)

zell <- zeller14[, zeller14$study_condition != "adenoma"]
nr <- nrow(zell)
exprs <- assay(zell, i = 1L)
datas <- data.frame(
    condition = lefser:::.numeric01(zell$study_condition),
    age = lefser:::.numeric01(zell$age_category),
    subjectID = zell$subjectID,
    t(exprs)
)

text <- mapply(function(x, y) {
    paste0(x, "\t", y)
    }, x = names(datas), y = sapply(datas, paste0, collapse = "\t")
)

## file for http://huttenhower.sph.harvard.edu/galaxy/ LEfSe
# writeLines(text, con = file("~/data/lefser/zelle14_150.txt"))
writeLines(text,
    con = file(paste0("~/data/lefser/zelle14_", nr, ".txt")))

## compared to results from lefser
loc <- Sys.getlocale("LC_COLLATE")
Sys.setlocale("LC_COLLATE", "C")
on.exit(Sys.setlocal("LC_COLLATE", loc))

# seed from
# https://github.com/SegataLab/lefse/blob/a4b3140c34b3abe5579a916a20c2a8659c3ac53c/lefse.py#L9
set.seed(1982)
res <- lefser(zell, groupCol = "study_condition", blockCol = "age_category")
lefserPlot(res)
