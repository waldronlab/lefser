# lefser 1.12.0

## Significant user-visible changes

* The `checkAbundances` argument in `lefser()` checks that data are as relative
abundances and warns if otherwise (@LiNk-NY @sdgamboa, #28)
* `relativeAb` helper function available to convert data (@LiNk-NY)
* Deprecate the `expr` argument and use `relab` (short for relative abundances)
* Add group labels to `lefserPlot` (@LiNk-NY #25, @asyakhl #31)
* 'Interoperating with `phyloseq`' section added to the vignette (#16)

# lefser 1.0.0

* `lefser` is an R/Bioconductor implementation of the `LEfSe` method for microbiome marker
discovery (https://doi.org/10.1186/gb-2011-12-6-r60)
* LEfSe uses the Kruskal-Wallis test, Wilcoxon-Rank Sum test, and Linear
Discriminant Analysis to find biomarkers in groups and sub-group blocks.
* `lefser` also implements the format of the LEfSe barplot of results
