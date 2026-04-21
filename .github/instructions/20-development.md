# Development Patterns

For complete Bioconductor and waldronlab standards, see:
- [Core Bioconductor standards](../../templates/bioconductor-development.md)
- [Waldronlab conventions](../../templates/waldronlab-standards.md)

## Package-Specific Patterns

### Function Organization

- `R/lefser.R`: Core LEfSe algorithm, mathematical modeling (`ldaFunction`, `filterKruskal`, `wilcox_pstats`), and the primary `lefser` function.
- `R/lefserClades.R`: Functions handling specific cladogram structuring, node manipulation, and taxonomy extraction from path strings.
- `R/lefserPlot*.R`: Visualization functions (`lefserPlot.R`, `lefserPlotClad.R`, `lefserPlotFeat.R`) that wrap `ggplot2` and `ggtree` to visualize LEfSe results.
- `R/utils.R`: Shared generic helper functions like `relativeAb`, `get_terminal_nodes`, `rowNames2RowData`, and plotting configuration (`.selectPalette`).
- Hidden helper functions are prefixed with a dot (e.g., `.trunc`, `.selectTaxRanks`, `.prepareDataHistogram`).

### Naming Conventions

- **Functions**: Generally use `camelCase` for primary and exported functions (e.g., `lefserClades`, `relativeAb`). Some generic wrappers use `snake_case` (e.g., `get_terminal_nodes`, `wilcox_pstats`). Internal helpers are prefixed with `.` (e.g., `.extractTips`).
- **Variables**: Descriptive, abbreviated names in `camelCase` or `snake_case` (e.g., `relab_sub`, `kruskal.threshold`).

## S4 Classes and Methods

While `lefser` does not define its own custom S4 classes, it relies heavily on core Bioconductor S4 classes as its primary data structures:
- `SummarizedExperiment::SummarizedExperiment`: The primary input structure containing abundance data (`relab`), sample metadata, and feature taxonomy.
- `S4Vectors::DataFrame` and `S4Vectors::SimpleList`: Used internally for manipulating metadata structures.
- Standard S4 methods actively used include `assay()`, `assays()`, `colData()`, and `rowData()`.

## Key Dependencies

SummarizedExperiment (Depends), coin, MASS, ggplot2, S4Vectors, dplyr, ggtree, mia, treeio (Imports)

## Code Style Notes

- Uses `roxygen2` for all function documentation.
- Extensively utilizes `ggplot2` and `ggtree` ecosystem for visualization outputs, manipulating aesthetics inside functional wrappers.
- Employs `tryCatch` and `withCallingHandlers` in the test suite to safely ignore specific expected upstream warnings (like `mia` or `MASS` deprecations and collinearity warnings).
