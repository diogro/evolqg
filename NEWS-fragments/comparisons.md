## Comparisons group (RandomSkewers, KrzCor, PCAsimilarity, MantelCor/MatrixCor, KrzProjection, MatrixDistance, DeltaZCorr, PrintMatrix)

### Changes in results

* `RandomSkewers()`, `KrzCor()`, `PCAsimilarity()` and `MatrixCor()` (and `MantelCor()` for samples) called with two `mcmc_sample` objects returned n^2 values (the whole first sample list compared against each sample of the second) or failed. They now compare sample i of the first with sample i of the second and return n values.
* `MatrixCor()` for an `mcmc_sample` against a single matrix now converts the samples (and the matrix) to correlations with `cov2cor`, as the sample-vs-sample method already did. Before, covariance samples were compared uncorrected.
* `MantelCor(landmark.dim = , withinLandmark = TRUE)` included the diagonal of the correlation matrices in the correlation (`lower.tri(x, diag = landmark.dim)`). The diagonal is now excluded; all other lower-triangle elements are used.
* `MantelCor()` with `landmark.dim` computed the p-value as `sum(null > observed)/permutations`, which could be 0. It is now `(sum(null >= observed) + 1)/(permutations + 1)`, the same as without landmarks.
* `MantelCor()` no longer prints vegan's "complete enumeration" messages; the p-values and the random stream of the non-landmark test are unchanged.
* `KrzCor()` and `KrzProjection()` with 2 or 3 traits: the default `ret.dim` (`round(n/2 - 1)`) was 0, giving `Inf` and an error. The default is now `max(1, round(n/2 - 1))`; it is unchanged for 4 or more traits.
* `RandomSkewers()` on a list: the random vectors are now drawn once per call and shared by all pairwise comparisons (and the null distribution for the probabilities is computed once). Results for a given seed differ from earlier versions; they are statistically equivalent. `RandomSkewers()` for a single pair of matrices keeps its random stream.
* `OverlapDist()` could return `NaN` for very similar matrices (`sqrt` of a tiny negative number); it now returns 0 in that case. `RiemannDist()` is computed from a symmetric problem (Cholesky whitening) and can no longer return complex values; results agree with the old ones to about 1e-12.
* `MatrixDistance()` (default method) now passes `...` (e.g. `iterations`) to the distance function.
* `DeltaZCorr()` on a list: with more than one skewer, `$probabilities` held the correlations of skewer 2 and skewers 3 and up were dropped. With no `cov.y` it now returns an n x n x k array (dimnames: matrix names, matrix names, skewer names; entry `[j, i, ]` compares matrices j and i, j > i). With one skewer the old list (`correlations`, `probabilities`) is kept. With `cov.y` it returns a data.frame with one column per skewer (plus `.id` for named lists) and no longer adds a column of `NA`.
* `DeltaZCorr()` accepts a vector as skewer and checks the dimensions of the skewers and matrices.
* `KrzProjection()` on lists with duplicated or missing names now makes the names unique (as `ComparisonMap` did) instead of returning a matrix of counts; `ComparisonMap()` accepts partially named lists and lists with a single matrix (1x1 result).
* Comparison methods for lists check that `cov.y` is a matrix with the same dimensions as the elements of the list.

### Bug fixes

* `MantelCor.mcmc_sample` passed `parallel` positionally into `...`; it is now passed by name.
* `PrintMatrix()` for lists no longer calls `dump()` and always closes the `sink()`, also on error.
* Exact comparisons of floating point diagonals (to decide if a matrix is a correlation matrix) now use a tolerance.

### Performance

* `KrzCor()` and `PCAsimilarity()` on lists compute the eigen decomposition once per matrix instead of once per pair (about 25x for 40 matrices of 10 traits).
* `MatrixCor()` on lists computes all correlations in one call to `cor()`.
* `KrzProjection()` (default) uses a single matrix product instead of nested `aaply` calls; the list method computes each eigen decomposition once and parallelises only the outer loop (about 100x for 40 matrices).
* `RandomSkewers()` on lists is about 7-10x faster (see above).
* `ComparisonMap()` / `SingleComparisonMap()` build their results directly instead of with `adply`/`acast`/`ldply`.
* The C++ functions `RS` and `delta_z_corr` take their arguments by const reference and check dimensions and `num_vectors`.

### Dependencies

* Code in these files no longer uses reshape2, vegan, mvtnorm or Matrix::bdiag (the Mantel test uses internal code with the same results and random stream).

### New arguments

* None. (`SingleComparisonMap()`, internal, gained `check.dim`.)
