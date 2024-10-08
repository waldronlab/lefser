
suppressPackageStartupMessages({
    library(MicrobiomeBenchmarkData)
    library(lefser)
})

## Example 1 - Strains
## Clades will be summarized at the species level and above
data("zeller14")
z14 <- zeller14[, zeller14$study_condition %in% c("control", "CRC")]
tn <- get_terminal_nodes(rownames(z14))
z14_tn <- z14[tn, ]
z14_tn_ra <- relativeAb(z14_tn)
z14_input <- rowNames2RowData(z14_tn_ra)

resCl <- lefserClades(relab = z14_input, classCol = "study_condition")
(ggt <- lefserPlotClad(resCl))

## Example 2 - OTUs
## Clades will be summarized at the genus level and above
tse <- getBenchmarkData("HMP_2012_16S_gingival_V35", dryrun = FALSE)[[1]]
tse_ra <- relativeAb(tse)
res_tse_cld <- lefserClades(relab = tse_ra, classCol = "body_subsite")

(ggt2 <- lefserPlotClad(res_tse_cld))

sessioninfo::session_info()
