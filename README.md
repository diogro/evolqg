EvolQG Package
======================

[![CRAN status](https://www.r-pkg.org/badges/version/evolqg)](https://CRAN.R-project.org/package=evolqg)
[![R-CMD-check](https://github.com/lem-usp/evolqg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lem-usp/evolqg/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/18354/lem-usp/EvolQG.svg)](https://zenodo.org/badge/latestdoi/18354/lem-usp/EvolQG)


Package for evolutionary quantitative genetics, described in Melo D, Garcia G, Hubbe A et al. EvolQG - An R package for evolutionary quantitative genetics [version 3; referees: 2 approved, 1 approved with reservations] F1000Research 2016, 4:925 (doi: [10.12688/f1000research.7082.3](https://doi.org/10.12688/f1000research.7082.3)).


Installation
============

From CRAN
---------

```R
install.packages("evolqg")
```

From GitHub
-----------

On Windows you first need [Rtools](https://cran.r-project.org/bin/windows/Rtools/), because the package has compiled code. Then:

```R
# install.packages("remotes")
remotes::install_github("lem-usp/evolqg")
```

Parallel computation
--------------------

Functions with a `parallel` argument run through whatever [foreach](https://cran.r-project.org/package=foreach) backend you register, for example:

```R
library(doParallel)
registerDoParallel(cores = 2)
MeanMatrixStatistics(RandomMatrix(10), parallel = TRUE)
```

Troubleshooting
---------------

All dependencies are available from CRAN. Please open an [issue](https://github.com/lem-usp/evolqg/issues) if the install fails for some reason.
