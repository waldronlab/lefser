# Instructions Index

Quick reference for AI assistants working with lefser.

## Quick Links

- [Overview](00-overview.md) - Package classification and key functions
- [Development](20-development.md) - Coding standards
- [Testing & Docs](30-testing-and-docs.md) - Testing and documentation
- [Vignettes](40-vignettes.md) - Vignette guide

## Package Metadata

- **Name**: lefser
- **Type**: Analysis Package
- **Version**: 1.21.7
- **Repository**: https://github.com/waldronlab/lefser

## Key Functions Quick Reference

**Data Processing Functions**: `lefser()`, `lefserClades()`, `relativeAb()`, `rowNames2RowData()`, `get_terminal_nodes()`
**Visualization Functions**: `lefserPlot()`, `lefserPlotClad()`, `lefserPlotFeat()`

## External Resources

- **Publication (Bioinformatics 2024):** [Lefser: Implementation of metagenomic biomarker discovery tool, LEfSe, in R](https://doi.org/10.1093/bioinformatics/btae707)
- **PubMed Central:** [PMC11665633](https://pmc.ncbi.nlm.nih.gov/articles/PMC11665633/)
- **Original LEfSe (Python/Galaxy):** [Segata et al. 2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3218848/)

## Quick Commands

```r
BiocManager::install("lefser")
library(lefser)
?lefser
browseVignettes("lefser")
```
