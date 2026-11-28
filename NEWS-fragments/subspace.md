## Changes in results

- `KrzSubspaceBootstrap()`: the null distribution is now a true permutation. In each replicate the pooled residuals are shuffled once and split into groups of the original sizes, so each residual is used exactly once. Before, each group drew its residuals independently from the pooled set, so one residual could land in several groups. The null P-matrices are now computed with the same residual degrees of freedom as the observed models (before they used `lm(~1)`, with a different df). `random` therefore changes; `observed` is unchanged for a given seed. The random stream also changes (all permutations are drawn first, in the master process, so parallel and serial runs use the same permutations).
- `LocalShapeVariables()` / `Center2MeanJacobianFast()`: the centring is now the Woods (2003) one, `logdet(J_i) - mean_i(logdet(J_i))`, which is the log determinant relative to the Karcher (log-Euclidean) mean. Before, `MeanMatrix()` used `chol()`, which reads only the upper triangle of the (non-symmetric) TPS jacobians, so the output equalled `logdet(J_i) - mean_i(logdet(sym_upper(J_i)))`. Each local variable shifts by a constant per volume; its variance and the covariance between volumes are unchanged (to about 1e-11) and the means are now exactly zero.
- `JacobianArray()`: the jacobian of the TPS is now computed in closed form instead of with `numDeriv::jacobian()`. Values change at about 1e-10 (the numerical error of the old method) and it is about 500 times faster. Arguments passed in `...` (for numDeriv) are ignored with a warning.
- `KrzSubspace()`: with p = 2 or 3 traits the default `k` is now 1 (it was 0, which failed). `k_angles` for identical matrices is 0 instead of NaN.

## Bug fixes

- `KrzSubspace()`: `k = 1` no longer errors; `k > p` gives a clear error; angles no longer return NaN from rounding above 1. The documented output name is `k_eVec_H`.
- `KrzSubspaceBootstrap()`: the default is `k = 3`, documented. `KrzSubspaceDataFrame(n = 1)` no longer errors; the documentation now says `mean` is the midpoint of the HPD interval.
- `RSProjection()` and `PlotRSprojection()` work with two populations; a single population gives a clear error. The example uses `ncols`.
- `SRD()` works with fewer than 40 iterations. `SRD()` on a list with one matrix gives a clear error and a non-NULL `cov.y` with a list gives a warning. Only the pairwise comparisons run in parallel when `parallel = TRUE`.
- `PCScoreCorrelation()` checks that `length(taxons)` equals the number of means.
- `MultiMahalanobis()` works with a single mean (returns a 1 x 1 zero matrix).
- `LocalShapeVariables()` computes `cs` when only `cs` is missing and gives clear messages when `gpa` is missing.

## Performance

- `MultiMahalanobis()` is vectorised (about 1000 times faster, same values).
- `SRD()` iterations are vectorised (about 150 times faster, same values and random stream).
- `RSProjection()` no longer builds a num.vectors x num.vectors product per matrix (about 75 times faster).
- `KrzSubspaceBootstrap()` computes only the shared matrix H per sample (about 40 times faster); the example now runs in seconds.
- `TPS()` builds the kernel matrix without per-element loops.

## Dependencies

- Morpho is now only suggested: `LocalShapeVariables()` needs it only when both `gpa` and `cs` are NULL. `cSize` and `arrMean3` are computed in base R.
- numDeriv, expm (in these functions), coda and reshape2 are no longer used by `KrzSubspace*`, `RSProjection`, `PlotRSprojection`, `JacobianArray`, `Center2MeanJacobianFast`, `MultiMahalanobis`.
- `aes_string()` replaced by `.data` in `PCScoreCorrelation()` and `PlotRSprojection()`.

## New arguments

- `LocalShapeVariables(verbose = TRUE)`: progress is reported with `message()` (it used `cat`/`print` and an always-on progress bar).
- `MultiMahalanobis(parallel)` is kept but ignored.
