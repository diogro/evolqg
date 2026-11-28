# evolqg 1.0-0

This release fixes a long list of bugs found in a review of the whole package,
makes many functions much faster, and cuts the dependencies roughly in half.

## Read this first

* **Some results change.** Every function whose output differs from 0.4-3 for
  the same input and seed is listed under "Changes in results" below, with the
  reason. Most are bug fixes. Functions not listed there give the same results
  as before (checked against 0.4-3 with fixed seeds).
* **Fewer dependencies.** evolqg no longer imports vegan, coda, MCMCpack,
  igraph, reshape2, expm, numDeriv, mvtnorm or Morpho. Morpho is now only
  suggested, for `LocalShapeVariables()` when no Procrustes alignment is
  supplied. The functions taken from these packages are reimplemented
  internally and give the same results, including the same random numbers for
  a given seed (Mantel tests, inverse Wishart draws, multivariate normal draws,
  HPD intervals). A fresh install now needs 25 other packages instead of 79,
  and no longer pulls in rgl and its system libraries.
* **plyr is still attached** by `library(evolqg)`, so scripts that use plyr
  functions without loading it keep working. `parallel = TRUE` still runs
  through whatever foreach backend you register (e.g. `doParallel`).
* **`LModularity()` is defunct** (see `?"evolqg-defunct"`).
* Tests use testthat 3rd edition; repository moved to
  <https://github.com/diogro/evolqg>.

## Changes in results

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
* `Rarefaction()` now passes `replace` to `RarefactionStat()`. Before, `replace = TRUE` was silently ignored and samples were always drawn without replacement. Results for `replace = TRUE` change.
* `CalcRepeatability()` estimates the among-individual variance with the Lessells & Boag (1987) coefficient n0 = (N - sum(n_i^2)/N)/(a - 1) instead of a hardcoded 2. For two measurements per individual n0 = 2 and results are unchanged; for unbalanced designs or more than two measurements the repeatabilities change. A single-column or vector `ind.data` is now accepted.
* `RarefactionStat()` and `Rarefaction()` now pass `...` to `ComparisonFunc`, as documented (previously dropped).
* `CalcEigenVar(sd = TRUE, sample = n)` returns `NA` with a warning, instead of `NaN`, when the observed variance is below its expected value. With `keep.positive = TRUE` the number of eigenvalues used throughout is the number of positive eigenvalues (now documented; unchanged).
* `BootstrapStat()`, `BootstrapRep()` and `MonteCarloStat()` with `iterations = 0` return zero rows (they returned 2 rows from `1:0`).
* `EigenTensorDecomposition()`: results differ from 0.4. The covariance tensor Sigma was built from the raw matrix elements, while the projections were taken in the log-centred space (`logm(M^-1/2 A M^-1/2)`, M the Riemannian mean), so eigentensors, eigenvalues and projections did not describe the same thing. Everything now lives in the log-centred (affine-invariant) space: the eigentensors are orthonormal in that space, `values` are the variances of the `projection` columns, and the projections are uncorrelated, mean-zero PCA-like scores. `ProjectMatrix()` and `RevertMatrix()` use the same space and are exact inverses of each other (with `scaled = FALSE`). Output names are unchanged; `mean.sqrt` now keeps the trait dimnames. With only two matrices the function no longer fails ("Invalid margin").
* `MeanMatrix()`: two typos relative to the Bini and Iannazzo (2013) algorithm were fixed (`beta` summed `ch` instead of `dh`, and `T + t(T)/2` instead of `(T + t(T))/2`). Results change at the level of the tolerance (about 1e-10), and ill-conditioned inputs that took minutes now converge in milliseconds. When the iteration gradient grows the previous (better) iterate is returned; a single matrix (m = 1) is returned as is; unnamed input no longer gets spurious dimnames.
* `RandomMatrix(LKJ = TRUE)` / LKJ correlations: the old generator capped first-row correlations at |r| <= 0.5 and used Beta(shape, shape) at every level, which is not an LKJ distribution. It now uses the vine method (Lewandowski et al. 2009) with `shape` as the LKJ parameter eta, so each correlation has a Beta(eta - 1 + n/2, eta - 1 + n/2) marginal on (-1, 1). Seeded LKJ output changes. The default partial-correlation method keeps exactly the same random stream and output.
* `CalculateMatrix()`: residuals are no longer re-centred, so models without an intercept give the correct matrix (`crossprod(residuals) / df.residual`). Identical (to 1e-16) for models with an intercept. A univariate `lm` returns a 1 x 1 matrix instead of `numeric(0)`.
* `RelativeEigenanalysis()`: the eigenvalues are always real (the old `eigen(cov.x %*% solve(cov.y))` returned complex values for many pairs of 30-trait matrices). The eigenvalues are computed from a Cholesky whitening with symmetric `eigen()`, and the eigenvectors are back-transformed; they keep the old scaling (unit length, same direction as the eigenvectors of `cov.x %*% solve(cov.y)`), with arbitrary sign.
* `MantelModTest()` / `TestModularity()`: with `landmark.dim` supplied, `withinLandmark` is now honoured. It used to be hardcoded to `FALSE`, so the documented option did nothing. The check that hypotheses are binary now uses `!all()` instead of `!any()`, so matrices mixing 0/1 with other values are rejected. Extra arguments in `...` are now passed on by the list method. Hypotheses with no within-module pairs (all zeros, single-trait modules) or with all trait pairs in modules now give a warning, since their correlation or AVG ratio is undefined.
* `MantelModTest()` / `TestModularity()` (list method, no `landmark.dim`): all hypotheses are now tested with one shared set of permutations, which is about 9 times faster. `rho` and the AVG columns are unchanged. The p-values (`Probability`) change under a fixed seed because the random stream is different and hypotheses now share permutations. A single hypothesis passed to `MantelModTest()` gives the same p-value as before for the same seed.
* `MINT()`: the best partition is now recovered from a stored column mask instead of splitting the hypothesis name on "_". It no longer fails when hypothesis names contain "_" or when the "null" hypothesis ranks first (it then warns and returns an identity hypothesis matrix). The `sig` column (fraction of Monte Carlo gammas below the observed gamma; values near 1 mean a poor fit) is unchanged. A new column `p.value` is added when `significance = TRUE`: (number of simulated gammas >= observed + 1)/(`iterations` + 1), with simulated and observed matrices compared on the correlation scale. High values mean the data are consistent with the hypothesis. Both columns are documented.
* `JackKnifeMINT()`: the upper limits (`upper.raw`, `upper.corrected`) used the 9.75% quantile (`0.0975`) instead of the 97.5% quantile. They are now correct.
* `PhyloMantel()`: the permuted matrix is now matched to the tips by name, so matrices in a different order from `tree$tip.label` give correct results (matrices are reordered to the tree, and missing names are an error). The null distribution now uses the supplied `ComparisonFunc` (it always used `cor` of the lower triangles). Triangular input (such as the output of `MatrixCor()` on a list) is symmetrised before permuting. The p-value is now (number of null values >= observed + 1)/(`permutations` + 1) instead of count/`permutations`.
* `MultivDriftTest()`: draws are now vectorised, which is about 400 times faster. The random stream changes, so results under a fixed seed differ slightly. It no longer fails with a single trait.
* `DriftTest()`: non-positive eigenvalues of the within-group matrix used to give NaN that was dropped silently from the regression. A warning is now given.
* `TreeDriftTest()`, `PlotTreeDriftTest()`: documentation now says "at least 4" descendant tips, which is what the code always required.
* `CreateHypotMatrix()`, `Partition2HypotMatrix()`: `Partition2HypotMatrix()` always returns a matrix (it returned a vector for a single trait or a single group) and errors on `NA`. `MINT()` and `JackKnifeMINT()` error above 20 hypotheses (the 32-bit mask used to break silently above 31) and warn above 15, because 2^n combinations are tested.
* `LModularity()` was removed (see Dependencies).
* `KrzSubspaceBootstrap()`: the null distribution is now a true permutation. In each replicate the pooled residuals are shuffled once and split into groups of the original sizes, so each residual is used exactly once. Before, each group drew its residuals independently from the pooled set, so one residual could land in several groups. The null P-matrices are now computed with the same residual degrees of freedom as the observed models (before they used `lm(~1)`, with a different df). `random` therefore changes; `observed` is unchanged for a given seed. The random stream also changes (all permutations are drawn first, in the master process, so parallel and serial runs use the same permutations).
* `LocalShapeVariables()` / `Center2MeanJacobianFast()`: the centring is now the Woods (2003) one, `logdet(J_i) - mean_i(logdet(J_i))`, which is the log determinant relative to the Karcher (log-Euclidean) mean. Before, `MeanMatrix()` used `chol()`, which reads only the upper triangle of the (non-symmetric) TPS jacobians, so the output equalled `logdet(J_i) - mean_i(logdet(sym_upper(J_i)))`. Each local variable shifts by a constant per volume; its variance and the covariance between volumes are unchanged (to about 1e-11) and the means are now exactly zero.
* `JacobianArray()`: the jacobian of the TPS is now computed in closed form instead of with `numDeriv::jacobian()`. Values change at about 1e-10 (the numerical error of the old method) and it is about 500 times faster. Arguments passed in `...` (for numDeriv) are ignored with a warning.
* `KrzSubspace()`: with p = 2 or 3 traits the default `k` is now 1 (it was 0, which failed). `k_angles` for identical matrices is 0 instead of NaN.

## Bug fixes

* `MantelCor.mcmc_sample` passed `parallel` positionally into `...`; it is now passed by name.
* `PrintMatrix()` for lists no longer calls `dump()` and always closes the `sink()`, also on error.
* Exact comparisons of floating point diagonals (to decide if a matrix is a correlation matrix) now use a tolerance.
* `BootstrapStat()`, `BootstrapRep()`, `BootstrapR2()`, `MonteCarloStat()`, `MonteCarloRep()` and `MonteCarloR2()` no longer loop forever, hiding the error, when the comparison errors or returns NA (for example one NA in the data, or a bad extra argument). They stop up front if the observed statistic has NAs, redraw a failing sample up to 100 times (warning once if any were redrawn) and then stop with the last error message.
* `RarefactionStat()` and `Rarefaction()` no longer loop forever when `StatFunc` always warns (a constant column with `correlation = TRUE`, or a deprecated function). The observed statistic is checked first (NA is an error; a warning there disables redraws), samples that warn are redrawn up to 100 times and then kept with one warning. Rarefaction needs at least 3 individuals (a clear error instead of a counting-down `3:n`); the documentation now says sizes start at 3, not 2.
* `parallel = TRUE` now gives the same populations as the serial run: bootstrap indices, Monte Carlo populations and rarefaction samples are all drawn in the main process in the serial order, and only the statistics are computed in parallel. Previously the draws happened in the workers and were not reproducible. Comparison functions that use random numbers themselves (e.g. `RandomSkewers`) still differ between workers.
* `CalcR2CvCorrected()` for `lm` models with more than one predictor failed ("object 'fac.s' not found"). Individuals are now grouped by the interaction of all factors. Its `parallel` argument, silently ignored before, now works; bootstrap indices are drawn in the main process. A resample with a single-individual group is still redrawn, but now stops with an error after 100 attempts in a row.
* `AlphaRep()`, `CalcICV()` and `MonteCarloStat()` test for a correlation matrix with a tolerance instead of exact floating-point equality of the diagonal.
* `BootstrapRep()` gained `...`, passed to `ComparisonFunc` like `MonteCarloRep()`.
* Internal use of the deprecated `CalcR2()` removed: `BootstrapR2()`, `MonteCarloR2()` and `CalcR2CvCorrected()` no longer emit one deprecation warning per iteration. `CalcR2()` itself stays deprecated.
* `PlotRarefaction()` uses `aes(.data$...)` instead of the deprecated `aes_string()`.
* `Autonomy()` and `ConditionalEvolvability()`: when the matrix was not positive definite and the `nearPD` fallback was used, the Cholesky factor of the nearest matrix was used instead of the matrix itself, giving nonsense (even negative conditional evolvabilities). The nearest positive definite matrix is now used, and numerically singular matrices that pass `chol()` are also sent through `nearPD`. A warning is still issued.
* `MeanMatrixStatistics()`: the critical value used `iterations` instead of the number of columns of a user-supplied `beta.mat` (printing `NA`); fixed. Matrices with one trait no longer error.
* `BayesianCalculateMatrix()`: a univariate `lm` no longer errors.
* `ExtendMatrix()`: matrices with fewer than 7 traits now give an error (they used to return junk), and a warning is given when no eigenvalue is extended because no gradient variance falls below `var.cut.off`.
* `RandomMatrix(1)` and `RandomMatrix(2)` (partial-correlation method) no longer fail.
* `maximum.norm()` (internal) had a typo that made it fail for vectors.
* Documentation of `CalculateMatrix()`, `BayesianCalculateMatrix()` (the "MLE" is the unbiased residual estimator, MAP is S_N/(nu + N), which is not the mode of the sampled inverse Wishart) and `RandomMatrix(shape =)` corrected; the wiki reference was replaced by explanations in the help pages, which also point to the stanimal package for fully Bayesian models.
* `PhyloCompare()`: children of each node were looked up by node number in `node.data`, but `PhyloW()` returns a list named by node in a different order (tips in `tip.data` order, then internal nodes in postorder). On the iris example tree the root was compared with itself (0.673 where the correct value is 1.010). Children are now looked up by name. Polytomies used to compare only the first two children; trees that are not binary now give an error.
* `PhyloW()`: node labels that were missing, empty or repeated (e.g. support values) made ancestral matrices overwrite each other silently. Internal nodes are now named by node number, with a warning, in that case. A named `tip.sample.size` keeps its names and is matched to the tips (it used to be renamed with `names(tip.data)`).
* `TreeDriftTest()` / `PlotTreeDriftTest()` failed on trees with node labels. Nodes are now handled by number internally and named by label (or number).
* `PhyloMantel()`, `JackKnifeMINT()`, `MINT()`, `MantelModTest()`, `MultivDriftTest()`: see "Changes in results".
* The internal `gamma()` in MINT masked `base::gamma()` inside the package; it is now `.mint_gamma()` and computed with `sum((S - S0)^2)` (12 times faster).
* `KrzSubspace()`: `k = 1` no longer errors; `k > p` gives a clear error; angles no longer return NaN from rounding above 1. The documented output name is `k_eVec_H`.
* `KrzSubspaceBootstrap()`: the default is `k = 3`, documented. `KrzSubspaceDataFrame(n = 1)` no longer errors; the documentation now says `mean` is the midpoint of the HPD interval.
* `RSProjection()` and `PlotRSprojection()` work with two populations; a single population gives a clear error. The example uses `ncols`.
* `SRD()` works with fewer than 40 iterations. `SRD()` on a list with one matrix gives a clear error and a non-NULL `cov.y` with a list gives a warning. Only the pairwise comparisons run in parallel when `parallel = TRUE`.
* `PCScoreCorrelation()` checks that `length(taxons)` equals the number of means.
* `MultiMahalanobis()` works with a single mean (returns a 1 x 1 zero matrix).
* `LocalShapeVariables()` computes `cs` when only `cs` is missing and gives clear messages when `gpa` is missing.

## Performance

* Loops that were slow in R now run in C++ (`src/fast_utils.cpp`). Random
  numbers are still drawn through R's generators in the same order, so seeded
  results are unchanged:
  * `RandomMatrix()`: both correlation generators (partial correlations, the
    default, and `LKJ = TRUE`) are 5-10 times faster. Values agree with earlier
    versions to about 1e-15; the old R loop's last digits depended on the BLAS.
  * Mantel permutation tests (`MantelCor()`, with or without `landmark.dim`,
    `MantelModTest()` and `TestModularity()`) compute the correlation for every
    permutation in one C++ call, about 2-2.5 times faster. Permutations and
    p-values are the same.
  * `MINT()` and `JackKnifeMINT()` build and deduplicate the combinations of
    hypotheses in C++ (about 5 times faster with 13 hypotheses).
* `KrzCor()` and `PCAsimilarity()` on lists compute the eigen decomposition once per matrix instead of once per pair (about 25x for 40 matrices of 10 traits).
* `MatrixCor()` on lists computes all correlations in one call to `cor()`.
* `KrzProjection()` (default) uses a single matrix product instead of nested `aaply` calls; the list method computes each eigen decomposition once and parallelises only the outer loop (about 100x for 40 matrices).
* `RandomSkewers()` on lists is about 7-10x faster (see above).
* `ComparisonMap()` / `SingleComparisonMap()` build their results directly instead of with `adply`/`acast`/`ldply`.
* The C++ functions `RS` and `delta_z_corr` take their arguments by const reference and check dimensions and `num_vectors`.
* `MonteCarloStat()` computes the SVD of the covariance matrix once instead of at every iteration (about 2.5x faster for cheap statistics). Same random stream and bit-identical results.
* `BootstrapStat()` and `RarefactionStat()` convert the data to a matrix once (when `StatFunc` is `cov` or `cor`) and resample row indices with `lapply`, instead of `alply`/`ldply` on a data frame (about 3x faster). Same random stream and identical results; custom `StatFunc`s still receive data frames.
* `CalcR2CvCorrected()` computes group coefficients of variation with `rowsum` and uses `eigen(only.values = TRUE)`.
* `MeanMatrixStatistics()` and `Autonomy()`, `ConditionalEvolvability()`, `Evolvability()`, `Flexibility()`, `Constraints()`, `Respondability()`, `Pc1Percent()`: no n x n intermediates, a single Cholesky with `backsolve`, and no `Matrix` coercion; values agree with 0.4 to about 1e-12 and the random stream is unchanged. `parallel = TRUE` still uses the registered foreach backend.
* `MeanMatrix()`: base R symmetric eigendecompositions instead of plyr and `Matrix::Schur` (about 20x), and many fewer iterations thanks to the typo fixes.
* `EigenTensorDecomposition()`: tensor algebra instead of nested plyr calls.
* `DriftTest()` is about 20 times faster with `plot = FALSE`; `CombineHypot()` (used by `MINT()` and `JackKnifeMINT()`) hashes the lower triangles instead of comparing every pair (about 10 times faster at 10 hypotheses); `MantelModTest()` list method shares permutations (about 9 times faster); `MultivDriftTest()` is vectorised; `PhyloMantel()` computes the cophenetic distances once.
* `MultiMahalanobis()` is vectorised (about 1000 times faster, same values).
* `SRD()` iterations are vectorised (about 150 times faster, same values and random stream).
* `RSProjection()` no longer builds a num.vectors x num.vectors product per matrix (about 75 times faster).
* `KrzSubspaceBootstrap()` computes only the shared matrix H per sample (about 40 times faster); the example now runs in seconds.
* `TPS()` builds the kernel matrix without per-element loops.

## Dependencies

* Code in these files no longer uses reshape2, vegan, mvtnorm or Matrix::bdiag (the Mantel test uses internal code with the same results and random stream).
* `reshape2` and `mvtnorm` are no longer used by these functions (`PlotRarefaction`, `Rarefaction`, `MonteCarloStat`, `BootstrapStat`). The `plyr` calls in these files were replaced by base R; `plyr` is still attached through Depends.
* `expm` is no longer used in `EigenTensorDecomposition()`, `ProjectMatrix()`, `RevertMatrix()` and `RelativeEigenanalysis()` (internal symmetric matrix functions are used instead); `MCMCpack::riwish` is replaced by an internal draw with the same random stream; `Matrix::Schur`, `Matrix::Matrix`, `solve`, `chol` and `diag` are no longer imported (only `nearPD`).
* `LModularity()` was removed because modularity maximisation is unreliable. It is now defunct and gives an error pointing to `igraph::cluster_optimal()`, `TestModularity()` and `MINT()`. This removes the igraph dependency (documented in the `evolqg-defunct` help topic).
* vegan, mvtnorm and igraph are no longer used by the modularity, drift and phylogeny functions (`.mantel_test()` and `.rmvnorm()` replace them). `ape::` is used explicitly. `ggplot2::aes_string()` was replaced by `aes(.data$x)`; `.data` is imported from ggplot2.
* Morpho is now only suggested: `LocalShapeVariables()` needs it only when both `gpa` and `cs` are NULL. `cSize` and `arrMean3` are computed in base R.
* numDeriv, expm (in these functions), coda and reshape2 are no longer used by `KrzSubspace*`, `RSProjection`, `PlotRSprojection`, `JacobianArray`, `Center2MeanJacobianFast`, `MultiMahalanobis`.
* `aes_string()` replaced by `.data` in `PCScoreCorrelation()` and `PlotRSprojection()`.

## New arguments

* `AlphaRep(method = c("cheverud", "squared", "exact"))`, as the last argument. `"cheverud"` (default, unchanged result) uses `(1 - mean(r)^2)/(N - 2)`; `"squared"` uses `(1 - mean(r^2))/(N - 2)`; `"exact"` uses `mean((1 - r^2)^2)/(N - 1)`.
* `BootstrapRep(...)`: additional arguments passed to `ComparisonFunc`.
* `CalcR2CvCorrected()` for `lm` objects: `parallel`.
* `DriftTest()`: `sample.sizes` (from a contributor, pull request #4) corrects the between-group variances for sampling error of the means. When `NULL` (default) results are unchanged. `tol` drops eigenvalues below it from the regression (default `NULL`, none dropped), and `plot = FALSE` skips building the ggplot (default `TRUE`).
* `TreeDriftTest()`: `plot = FALSE` skips building the per-node plots (default `TRUE`).
* `JackKnifeMINT()`: `parallel` runs the jackknife samples through a registered foreach backend. `...` is now ignored (it used to be passed to `plyr::raply()`).
* `MantelModTest()`: `perm.matrix` takes a precomputed permutation matrix.
* `LocalShapeVariables(verbose = TRUE)`: progress is reported with `message()` (it used `cat`/`print` and an always-on progress bar).
* `MultiMahalanobis(parallel)` is kept but ignored.

# evolqg 0.4-3

* Use `\doi{}` markup in the `ratones` documentation.
