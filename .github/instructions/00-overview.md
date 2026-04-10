# lefser Overview

## Classification
- **Type**: Analysis Package
- **Version**: 1.21.7

## Purpose

lefser is the R implementation of the popular microbiome biomarker discovery tool, LEfSe. It uses the Kruskal-Wallis test, Wilcoxon-Rank Sum test, and Linear Discriminant Analysis to find biomarkers from two-level classes (and optional sub-classes).

## Key Functions

**Data Processing Functions** (5):
- `lefser()` - Main function for LEfSe analysis
- `lefserClades()` - Agglomerate features and perform LEfSe analysis across taxonomic ranks
- `relativeAb()` - Convert feature abundances to relative abundances
- `rowNames2RowData()` - Move rownames to RowData
- `get_terminal_nodes()` - Filter to terminal nodes (leaves)

**Visualization Functions** (3):
- `lefserPlot()` - Plot LEfSe results
- `lefserPlotClad()` - Plot cladogram of LEfSe results
- `lefserPlotFeat()` - Plot feature distributions for LEfSe results

## Quick Start

```r
library(lefser)

# Load example data
data("zeller14")
# Subsetting and filtering to terminal taxonomic nodes to avoid collinearity
z14 <- zeller14[, zeller14$study_condition != "adenoma"]
tn <- get_terminal_nodes(rownames(z14))
z14tn <- z14[tn, ]

# Calculate relative abundance
z14tn_ra <- relativeAb(z14tn)

# Run LEfSe analysis
res <- lefser(z14tn_ra, classCol = "study_condition")

# Plot results
lefserPlot(res)
```

## Key Concepts

- **LEfSe Algorithm**: Linear discriminant analysis Effect Size. A biomarker discovery algorithm that identifies features (like microbial taxa) that are significantly different between biological groups by estimating the effect size (biological consistency).
- **Statistical Pipeline**: Uses the Kruskal-Wallis test (for differential abundance between primary classes), Wilcoxon Rank-Sum test (optional, for consistency across subclasses), and Linear Discriminant Analysis (LDA, to estimate the effect size).
- **Cladogram**: A hierarchical, tree-like visualization that illustrates the taxonomic relationships of the identified biomarkers.
- **Collinearity**: When working with hierarchical taxonomic data, non-terminal nodes (e.g., genus) are linearly dependent on their subset terminal nodes (e.g., species), which can cause MASS::lda to fail. The `get_terminal_nodes` function solves this.

## Data Sources

The package provides a built-in dataset `zeller14` for testing and examples. This is derived from a colorectal cancer (CRC) metagenomic study (Zeller et al. 2014) containing microbiome abundance data.
