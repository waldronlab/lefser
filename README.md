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

## To Cite *lefser* software package

To cite `lefser` in publications, use:

> <p>
> Asya Khleborodova, Samuel D Gamboa-Tuz, Marcel Ramos, Nicola Segata, Levi Waldron, Sehyun Oh,
> <em>Lefser</em>: Implementation of metagenomic biomarker 
> discovery tool, <em>LEfSe</em.>, in R, <em>Bioinformatics</em>, 2024;, btae707, <a 
> href="https://doi.org/10.1093/bioinformatics/btae707">doi:10.1093/bioinformatics/btae707</a>,
> <a href="https://doi.org/10.1093/bioinformatics/btae707">https://doi.org/10.1093/bioinformatics/btae707</a>.
> </p>

A BibTeX entry for LaTeX users is: 

     @article{10.1093/bioinformatics/btae707,
      author = {Khleborodova, Asya and Gamboa-Tuz, Samuel D and Ramos, Marcel and Segata, Nicola and Waldron, Levi and Oh, Sehyun},
      title = {Lefser: Implementation of metagenomic biomarker discovery tool, LEfSe, in R},
      journal = {Bioinformatics},
      pages = {btae707},
      year = {2024},
      month = {11},
      issn = {1367-4811},
      doi = {10.1093/bioinformatics/btae707},
      url = {https://doi.org/10.1093/bioinformatics/btae707},
      eprint = {https://academic.oup.com/bioinformatics/advance-article-pdf/doi/10.1093/bioinformatics/btae707/60811200/btae707.pdf},
    }