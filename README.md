<!-- badges: start -->
![build](https://github.com/waldronlab/lefser/workflows/build/badge.svg)
[![Codecov test coverage](https://codecov.io/gh/waldronlab/lefser/branch/devel/graph/badge.svg)](https://codecov.io/gh/waldronlab/lefser?branch=devel)
<!-- badges: end -->

## *lefser*: Run *LEfSe* in R
*lefser* is the R implementation of the Python package, Linear discriminant 
analysis (LDA) Effect Size (*[LEfSe][]*). *LEfSe* is the most widely used Python 
package and Galaxy module for metagenomic biomarker discovery and 
visualization ([Segata et al. 2011][]). *LEfSe* utilizes standard 
statistical significance tests along with supplementary tests that incorporate 
biological consistency and the relevance of effects to identity the features 
(e.g., organisms, clades, OTU, genes, or functions) that are most likely to 
account for differences between the two sample classes of interest, referred as
‘classes’. While *LEfSe* is widely used and available in different platform 
such as Galaxy UI and Conda, there is no convenient way to incorporate it in 
R-based workflows. Thus, we re-implement *LEfSe* as an R/Bioconductor package, 
*lefser*. Following the *LEfSe*‘s algorithm including Kruskal-Wallis test, 
Wilcoxon-Rank Sum test, and Linear Discriminant Analysis, with some 
modifications, *lefser* successfully reproduces and improves the original 
statistical method and the associated plotting functionality.

[LEfSe]: https://huttenhower.sph.harvard.edu/galaxy/
[Segata et al. 2011]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3218848/
