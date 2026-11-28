## Changes in results (resampling group)

* `Rarefaction()` now passes `replace` to `RarefactionStat()`. Before, `replace = TRUE` was silently ignored and samples were always drawn without replacement. Results for `replace = TRUE` change.
* `CalcRepeatability()` estimates the among-individual variance with the Lessells & Boag (1987) coefficient n0 = (N - sum(n_i^2)/N)/(a - 1) instead of a hardcoded 2. For two measurements per individual n0 = 2 and results are unchanged; for unbalanced designs or more than two measurements the repeatabilities change. A single-column or vector `ind.data` is now accepted.
* `RarefactionStat()` and `Rarefaction()` now pass `...` to `ComparisonFunc`, as documented (previously dropped).
* `CalcEigenVar(sd = TRUE, sample = n)` returns `NA` with a warning, instead of `NaN`, when the observed variance is below its expected value. With `keep.positive = TRUE` the number of eigenvalues used throughout is the number of positive eigenvalues (now documented; unchanged).
* `BootstrapStat()`, `BootstrapRep()` and `MonteCarloStat()` with `iterations = 0` return zero rows (they returned 2 rows from `1:0`).

## Bug fixes

* `BootstrapStat()`, `BootstrapRep()`, `BootstrapR2()`, `MonteCarloStat()`, `MonteCarloRep()` and `MonteCarloR2()` no longer loop forever, hiding the error, when the comparison errors or returns NA (for example one NA in the data, or a bad extra argument). They stop up front if the observed statistic has NAs, redraw a failing sample up to 100 times (warning once if any were redrawn) and then stop with the last error message.
* `RarefactionStat()` and `Rarefaction()` no longer loop forever when `StatFunc` always warns (a constant column with `correlation = TRUE`, or a deprecated function). The observed statistic is checked first (NA is an error; a warning there disables redraws), samples that warn are redrawn up to 100 times and then kept with one warning. Rarefaction needs at least 3 individuals (a clear error instead of a counting-down `3:n`); the documentation now says sizes start at 3, not 2.
* `parallel = TRUE` now gives the same populations as the serial run: bootstrap indices, Monte Carlo populations and rarefaction samples are all drawn in the main process in the serial order, and only the statistics are computed in parallel. Previously the draws happened in the workers and were not reproducible. Comparison functions that use random numbers themselves (e.g. `RandomSkewers`) still differ between workers.
* `CalcR2CvCorrected()` for `lm` models with more than one predictor failed ("object 'fac.s' not found"). Individuals are now grouped by the interaction of all factors. Its `parallel` argument, silently ignored before, now works; bootstrap indices are drawn in the main process. A resample with a single-individual group is still redrawn, but now stops with an error after 100 attempts in a row.
* `AlphaRep()`, `CalcICV()` and `MonteCarloStat()` test for a correlation matrix with a tolerance instead of exact floating-point equality of the diagonal.
* `BootstrapRep()` gained `...`, passed to `ComparisonFunc` like `MonteCarloRep()`.
* Internal use of the deprecated `CalcR2()` removed: `BootstrapR2()`, `MonteCarloR2()` and `CalcR2CvCorrected()` no longer emit one deprecation warning per iteration. `CalcR2()` itself stays deprecated.
* `PlotRarefaction()` uses `aes(.data$...)` instead of the deprecated `aes_string()`.

## Performance

* `MonteCarloStat()` computes the SVD of the covariance matrix once instead of at every iteration (about 2.5x faster for cheap statistics). Same random stream and bit-identical results.
* `BootstrapStat()` and `RarefactionStat()` convert the data to a matrix once (when `StatFunc` is `cov` or `cor`) and resample row indices with `lapply`, instead of `alply`/`ldply` on a data frame (about 3x faster). Same random stream and identical results; custom `StatFunc`s still receive data frames.
* `CalcR2CvCorrected()` computes group coefficients of variation with `rowsum` and uses `eigen(only.values = TRUE)`.

## Dependencies

* `reshape2` and `mvtnorm` are no longer used by these functions (`PlotRarefaction`, `Rarefaction`, `MonteCarloStat`, `BootstrapStat`). The `plyr` calls in these files were replaced by base R; `plyr` is still attached through Depends.

## New arguments

* `AlphaRep(method = c("cheverud", "squared", "exact"))`, as the last argument. `"cheverud"` (default, unchanged result) uses `(1 - mean(r)^2)/(N - 2)`; `"squared"` uses `(1 - mean(r^2))/(N - 2)`; `"exact"` uses `mean((1 - r^2)^2)/(N - 1)`.
* `BootstrapRep(...)`: additional arguments passed to `ComparisonFunc`.
* `CalcR2CvCorrected()` for `lm` objects: `parallel`.
