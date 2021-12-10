<!-- badges: start -->
![build](https://github.com/waldronlab/lefser/workflows/build/badge.svg)
[![Codecov test coverage](https://codecov.io/gh/waldronlab/lefser/branch/devel/graph/badge.svg)](https://codecov.io/gh/waldronlab/lefser?branch=devel)
<!-- badges: end -->

# lefser: Run LEfSE in R

`lefser` is the R implementation of the _LEfSe_ method for microbiome biomarker discovery[1]. The original software is likely the most widely-used method for biomarker discovery and plotting in microbiome studies, with ~5,000 citations as of the end of 2020. However, the original software is implemented in Python as a command-line tool and Galaxy module, and is not straightforward to implement in R. The method involves Kruskal-Wallis test, Wilcoxon-Rank Sum test, and Linear Discriminant Analysis to find biomarkers in groups and sub-group blocks. `lefser` closely reproduces the original statistical method and the associated barplot of results.


[1] Segata N, Izard J, Waldron L, Gevers D, Miropolsky L, Garrett WS, Huttenhower C: Metagenomic biomarker discovery and explanation. Genome Biol. 2011, 12:R60.  https://doi.org/10.1186/gb-2011-12-6-r60
