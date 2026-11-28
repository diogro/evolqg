## Changes in results

- `MantelModTest()` / `TestModularity()`: with `landmark.dim` supplied, `withinLandmark` is now honoured. It used to be hardcoded to `FALSE`, so the documented option did nothing. The check that hypotheses are binary now uses `!all()` instead of `!any()`, so matrices mixing 0/1 with other values are rejected. Extra arguments in `...` are now passed on by the list method. Hypotheses with no within-module pairs (all zeros, single-trait modules) or with all trait pairs in modules now give a warning, since their correlation or AVG ratio is undefined.
- `MantelModTest()` / `TestModularity()` (list method, no `landmark.dim`): all hypotheses are now tested with one shared set of permutations, which is about 9 times faster. `rho` and the AVG columns are unchanged. The p-values (`Probability`) change under a fixed seed because the random stream is different and hypotheses now share permutations. A single hypothesis passed to `MantelModTest()` gives the same p-value as before for the same seed.
- `MINT()`: the best partition is now recovered from a stored column mask instead of splitting the hypothesis name on "_". It no longer fails when hypothesis names contain "_" or when the "null" hypothesis ranks first (it then warns and returns an identity hypothesis matrix). The `sig` column (fraction of Monte Carlo gammas below the observed gamma; values near 1 mean a poor fit) is unchanged. A new column `p.value` is added when `significance = TRUE`: (number of simulated gammas >= observed + 1)/(`iterations` + 1), with simulated and observed matrices compared on the correlation scale. High values mean the data are consistent with the hypothesis. Both columns are documented.
- `JackKnifeMINT()`: the upper limits (`upper.raw`, `upper.corrected`) used the 9.75% quantile (`0.0975`) instead of the 97.5% quantile. They are now correct.
- `PhyloMantel()`: the permuted matrix is now matched to the tips by name, so matrices in a different order from `tree$tip.label` give correct results (matrices are reordered to the tree, and missing names are an error). The null distribution now uses the supplied `ComparisonFunc` (it always used `cor` of the lower triangles). Triangular input (such as the output of `MatrixCor()` on a list) is symmetrised before permuting. The p-value is now (number of null values >= observed + 1)/(`permutations` + 1) instead of count/`permutations`.
- `MultivDriftTest()`: draws are now vectorised, which is about 400 times faster. The random stream changes, so results under a fixed seed differ slightly. It no longer fails with a single trait.
- `DriftTest()`: non-positive eigenvalues of the within-group matrix used to give NaN that was dropped silently from the regression. A warning is now given.
- `TreeDriftTest()`, `PlotTreeDriftTest()`: documentation now says "at least 4" descendant tips, which is what the code always required.
- `CreateHypotMatrix()`, `Partition2HypotMatrix()`: `Partition2HypotMatrix()` always returns a matrix (it returned a vector for a single trait or a single group) and errors on `NA`. `MINT()` and `JackKnifeMINT()` error above 20 hypotheses (the 32-bit mask used to break silently above 31) and warn above 15, because 2^n combinations are tested.
- `LModularity()` was removed (see Dependencies).

## Bug fixes

- `PhyloCompare()`: children of each node were looked up by node number in `node.data`, but `PhyloW()` returns a list named by node in a different order (tips in `tip.data` order, then internal nodes in postorder). On the iris example tree the root was compared with itself (0.673 where the correct value is 1.010). Children are now looked up by name. Polytomies used to compare only the first two children; trees that are not binary now give an error.
- `PhyloW()`: node labels that were missing, empty or repeated (e.g. support values) made ancestral matrices overwrite each other silently. Internal nodes are now named by node number, with a warning, in that case. A named `tip.sample.size` keeps its names and is matched to the tips (it used to be renamed with `names(tip.data)`).
- `TreeDriftTest()` / `PlotTreeDriftTest()` failed on trees with node labels. Nodes are now handled by number internally and named by label (or number).
- `PhyloMantel()`, `JackKnifeMINT()`, `MINT()`, `MantelModTest()`, `MultivDriftTest()`: see "Changes in results".
- The internal `gamma()` in MINT masked `base::gamma()` inside the package; it is now `.mint_gamma()` and computed with `sum((S - S0)^2)` (12 times faster).

## Performance

- `DriftTest()` is about 20 times faster with `plot = FALSE`; `CombineHypot()` (used by `MINT()` and `JackKnifeMINT()`) hashes the lower triangles instead of comparing every pair (about 10 times faster at 10 hypotheses); `MantelModTest()` list method shares permutations (about 9 times faster); `MultivDriftTest()` is vectorised; `PhyloMantel()` computes the cophenetic distances once.

## Dependencies

- `LModularity()` was removed because modularity maximisation is unreliable. It is now defunct and gives an error pointing to `igraph::cluster_optimal()`, `TestModularity()` and `MINT()`. This removes the igraph dependency (documented in the `evolqg-defunct` help topic).
- vegan, mvtnorm and igraph are no longer used by the modularity, drift and phylogeny functions (`.mantel_test()` and `.rmvnorm()` replace them). `ape::` is used explicitly. `ggplot2::aes_string()` was replaced by `aes(.data$x)`; `.data` is imported from ggplot2.

## New arguments

- `DriftTest()`: `sample.sizes` (from a contributor, pull request #4) corrects the between-group variances for sampling error of the means. When `NULL` (default) results are unchanged. `tol` drops eigenvalues below it from the regression (default `NULL`, none dropped), and `plot = FALSE` skips building the ggplot (default `TRUE`).
- `TreeDriftTest()`: `plot = FALSE` skips building the per-node plots (default `TRUE`).
- `JackKnifeMINT()`: `parallel` runs the jackknife samples through a registered foreach backend. `...` is now ignored (it used to be passed to `plyr::raply()`).
- `MantelModTest()`: `perm.matrix` takes a precomputed permutation matrix.
