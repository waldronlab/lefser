# Testing and Documentation

For complete standards, see:
- [Core Bioconductor standards](../../templates/bioconductor-development.md)
- [Waldronlab conventions](../../templates/waldronlab-standards.md)

## Package-Specific Testing

### Test Organization

Tests are located in `tests/testthat/` with 3 test files (`test-lefser.R`, `test-lefserPlotClad.R`, `test-lefserPlotFeat.R`).

### Test Data

- **Location**: `data/zeller14.rda`
- **File types**: RDA files
- **Purpose**: Provides a demo dataset for examples and testing.

### Remote Data Testing

Not applicable.

### Running Tests

devtools::test()

## Package-Specific Documentation Patterns

### Function Categories

- **Core Analysis (`lefser`)**: Functions dealing directly with the statistical pipeline (Kruskal-Wallis, Wilcoxon, LDA).
- **Visualization (`lefserPlot`, `lefserPlotClad`, `lefserPlotFeat`)**: Wrappers returning `ggplot` or `ggtree` objects.
- **Data Utility (`relativeAb`, `get_terminal_nodes`)**: Preparing data and handling collinearity.

### Common Parameters

- `relab`: A `SummarizedExperiment` object with relative abundances.
- `classCol`: The name of the `colData` column defining the primary grouping variable.
- `subclassCol`: (Optional) The name of the `colData` column defining secondary/nested grouping variables.
- `kruskal.threshold`, `wilcox.threshold`, `lda.threshold`: Numeric thresholds corresponding to significance cutoffs.

## Common Testing Patterns

- **Environment Isolation**: Loading mock or example data into separate, temporary environments (e.g., using `new.env()`) for strict baseline testing.
- **Tolerance Checking**: Verifying statistical accuracy using `expect_equal(..., tolerance = tol)` since LDA scores compute continuous floating points.
- **Warning Suppression**: Wrapping known deprecation warnings inside `withCallingHandlers` coupled with `invokeRestart("muffleWarning")` and using `expect_no_warning()` strictly to prevent future un-handled regressions.
    - **TODO**: Remove the `withCallingHandlers` suppression for `mia`'s `sumCountsAcrossFeatures` in `test-lefserPlotClad.R` once `mia::splitByRanks` updates its internal implementation to `aggregateAcrossFeatures` and the upstream deprecation warning is resolved.
- **Cross-Platform Numerical Stability**: Mathematical dependencies (like `MASS::lda`) evaluate singular value decompositions using underlying system libraries (BLAS/LAPACK). This means identical collinear test data might evaluate to Rank 1 on Linux (due to float noise) but Rank 0 on macOS (throwing an error). Test data structures are deliberately augmented to break perfect mathematical symmetry to ensure cross-platform test stability without requiring `tryCatch` error absorption.
