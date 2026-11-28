## Changes in results

- `EigenTensorDecomposition()`: results differ from 0.4. The covariance tensor Sigma was built from the raw matrix elements, while the projections were taken in the log-centred space (`logm(M^-1/2 A M^-1/2)`, M the Riemannian mean), so eigentensors, eigenvalues and projections did not describe the same thing. Everything now lives in the log-centred (affine-invariant) space: the eigentensors are orthonormal in that space, `values` are the variances of the `projection` columns, and the projections are uncorrelated, mean-zero PCA-like scores. `ProjectMatrix()` and `RevertMatrix()` use the same space and are exact inverses of each other (with `scaled = FALSE`). Output names are unchanged; `mean.sqrt` now keeps the trait dimnames. With only two matrices the function no longer fails ("Invalid margin").
- `MeanMatrix()`: two typos relative to the Bini and Iannazzo (2013) algorithm were fixed (`beta` summed `ch` instead of `dh`, and `T + t(T)/2` instead of `(T + t(T))/2`). Results change at the level of the tolerance (about 1e-10), and ill-conditioned inputs that took minutes now converge in milliseconds. When the iteration gradient grows the previous (better) iterate is returned; a single matrix (m = 1) is returned as is; unnamed input no longer gets spurious dimnames.
- `RandomMatrix(LKJ = TRUE)` / LKJ correlations: the old generator capped first-row correlations at |r| <= 0.5 and used Beta(shape, shape) at every level, which is not an LKJ distribution. It now uses the vine method (Lewandowski et al. 2009) with `shape` as the LKJ parameter eta, so each correlation has a Beta(eta - 1 + n/2, eta - 1 + n/2) marginal on (-1, 1). Seeded LKJ output changes. The default partial-correlation method keeps exactly the same random stream and output.
- `CalculateMatrix()`: residuals are no longer re-centred, so models without an intercept give the correct matrix (`crossprod(residuals) / df.residual`). Identical (to 1e-16) for models with an intercept. A univariate `lm` returns a 1 x 1 matrix instead of `numeric(0)`.
- `RelativeEigenanalysis()`: the eigenvalues are always real (the old `eigen(cov.x %*% solve(cov.y))` returned complex values for many pairs of 30-trait matrices). The eigenvalues are computed from a Cholesky whitening with symmetric `eigen()`, and the eigenvectors are back-transformed; they keep the old scaling (unit length, same direction as the eigenvectors of `cov.x %*% solve(cov.y)`), with arbitrary sign.

## Bug fixes

- `Autonomy()` and `ConditionalEvolvability()`: when the matrix was not positive definite and the `nearPD` fallback was used, the Cholesky factor of the nearest matrix was used instead of the matrix itself, giving nonsense (even negative conditional evolvabilities). The nearest positive definite matrix is now used, and numerically singular matrices that pass `chol()` are also sent through `nearPD`. A warning is still issued.
- `MeanMatrixStatistics()`: the critical value used `iterations` instead of the number of columns of a user-supplied `beta.mat` (printing `NA`); fixed. Matrices with one trait no longer error.
- `BayesianCalculateMatrix()`: a univariate `lm` no longer errors.
- `ExtendMatrix()`: matrices with fewer than 7 traits now give an error (they used to return junk), and a warning is given when no eigenvalue is extended because no gradient variance falls below `var.cut.off`.
- `RandomMatrix(1)` and `RandomMatrix(2)` (partial-correlation method) no longer fail.
- `maximum.norm()` (internal) had a typo that made it fail for vectors.
- Documentation of `CalculateMatrix()`, `BayesianCalculateMatrix()` (the "MLE" is the unbiased residual estimator, MAP is S_N/(nu + N), which is not the mode of the sampled inverse Wishart) and `RandomMatrix(shape =)` corrected; the wiki reference was replaced by explanations in the help pages, which also point to the stanimal package for fully Bayesian models.

## Performance

- `MeanMatrixStatistics()` and `Autonomy()`, `ConditionalEvolvability()`, `Evolvability()`, `Flexibility()`, `Constraints()`, `Respondability()`, `Pc1Percent()`: no n x n intermediates, a single Cholesky with `backsolve`, and no `Matrix` coercion; values agree with 0.4 to about 1e-12 and the random stream is unchanged. `parallel = TRUE` still uses the registered foreach backend.
- `MeanMatrix()`: base R symmetric eigendecompositions instead of plyr and `Matrix::Schur` (about 20x), and many fewer iterations thanks to the typo fixes.
- `EigenTensorDecomposition()`: tensor algebra instead of nested plyr calls.

## Dependencies

- `expm` is no longer used in `EigenTensorDecomposition()`, `ProjectMatrix()`, `RevertMatrix()` and `RelativeEigenanalysis()` (internal symmetric matrix functions are used instead); `MCMCpack::riwish` is replaced by an internal draw with the same random stream; `Matrix::Schur`, `Matrix::Matrix`, `solve`, `chol` and `diag` are no longer imported (only `nearPD`).
