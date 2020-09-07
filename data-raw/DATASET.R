## code to prepare example dataset goes here
## can change dataset by replacing ZellerG_2014 with any other dataset

library(curatedMetagenomicData)
exampledata <-
  curatedMetagenomicData("ZellerG_2014.metaphlan_bugs_list.stool",
                         counts = TRUE,
                         dryrun = FALSE)[[1]]
exampledata <- as(exampledata, "SummarizedExperiment")
usethis::use_data(exampledata, overwrite = TRUE)
