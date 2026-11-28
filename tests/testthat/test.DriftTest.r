# The lm object keeps the function frame as its formula environment, so
# expect_equal() on whole results compares everything defined in the call
# (including the order of the means). Compare the contents instead.
.drift_content <- function(x) c(list(coefficients = coef(x$regression)), x[c(2:5)])

test_that("DriftTest returns resonable results",
{
  local_seed(1)
  means <- array(rnorm(40*10), c(10, 40)) 
  means.list <- alply(means, 1)
  cov.matrix <- RandomMatrix(40, 1, 1, 10, LKJ = FALSE)
  test.array <- DriftTest(means, cov.matrix, FALSE)
  test.list <- DriftTest(means.list, cov.matrix, FALSE)
  expect_equal(.drift_content(test.array), .drift_content(test.list))
  expect_type(test.array, "list")
  expect_s3_class(test.array[[1]], "lm")
  expect_type(test.array[[2]], "double")
  expect_true(is.matrix(test.array[[2]]))
  expect_type(test.array[[3]], "double")
  expect_type(test.array[[4]], "double")
  expect_type(test.array[[5]], "logical")
  expect_s3_class(test.array[[6]], c("gg", "ggplot"))
  expect_equal(names(test.array), c("regression",
                                    "coefficient_CI_95",
                                    "log.between_group_variance",
                                    "log.W_eVals",
                                    "drift_rejected",
                                    "plot"))
  expect_equal(length(test.array$log.between_group_variance), dim(means)[2])
  expect_equal(length(test.array$log.W_eVals), dim(means)[2])
  expect_equal(test.array$log.W_eVals, log(eigen(cov.matrix)$values))
})

test_that("DriftTest plot = FALSE skips the plot only", {
  local_seed(1)
  means <- array(rnorm(20*10), c(10, 20))
  cov.matrix <- RandomMatrix(20, 1, 1, 10, LKJ = FALSE)
  with.plot <- DriftTest(means, cov.matrix, FALSE)
  no.plot <- DriftTest(means, cov.matrix, FALSE, plot = FALSE)
  expect_equal(names(no.plot), names(with.plot))
  expect_null(no.plot$plot)
  expect_equal(.drift_content(no.plot), .drift_content(with.plot))
  # show.plot = TRUE still draws
  pdf(NULL); on.exit(dev.off())
  expect_s3_class(DriftTest(means, cov.matrix, TRUE, plot = FALSE)$plot, "ggplot")
})

test_that("DriftTest warns about non-positive-definite matrices and has a tol option", {
  local_seed(2)
  p <- 8
  means <- array(rnorm(12*p), c(12, p))
  W <- RandomMatrix(p, 1, 1, 10, LKJ = FALSE)
  ev <- eigen(W)
  bad <- ev$vectors %*% diag(c(ev$values[1:6], 1e-9, -0.5)) %*% t(ev$vectors)
  bad <- (bad + t(bad))/2
  expect_warning(res <- DriftTest(means, bad, FALSE, plot = FALSE), "not positive definite")
  expect_equal(sum(is.na(res$log.W_eVals)), 1)
  expect_equal(length(res$regression$residuals), 7)
  # tol drops the near-zero and negative eigenvalues, with no warning
  expect_warning(res.tol <- DriftTest(means, bad, FALSE, tol = 1e-6, plot = FALSE), NA)
  expect_equal(length(res.tol$regression$residuals), 6)
  psd <- ev$vectors %*% diag(c(ev$values[1:7], 1e-12)) %*% t(ev$vectors)
  psd <- (psd + t(psd))/2
  expect_warning(res.tol <- DriftTest(means, psd, FALSE, tol = 1e-6, plot = FALSE), NA)
  expect_equal(length(res.tol$regression$residuals), 7)
  expect_equal(length(res.tol$log.W_eVals), p)
  # default behaviour is unchanged for positive definite matrices
  expect_warning(DriftTest(means, W, FALSE, plot = FALSE), NA)
})

.drift_tree_data <- function(tree) {
  mean.list <- lapply(tree$tip.label, function(x) rnorm(5))
  names(mean.list) <- tree$tip.label
  cov.matrix.list <- RandomMatrix(5, length(tree$tip.label))
  names(cov.matrix.list) <- tree$tip.label
  list(means = mean.list, covs = cov.matrix.list)
}

test_that("TreeDriftTest returns resonable results",
{
  local_seed(3)
  library(ape)
  data(bird.orders)
  tree <- bird.orders
  d <- .drift_tree_data(tree)
  w.cov <- PhyloW(tree, d$covs)$'24'
  test.list <- TreeDriftTest(tree, d$means, d$covs)
  expect_type(test.list, "list")
  expect_equal(length(test.list), 13)
  expect_equal(.drift_content(test.list[[length(test.list)]]),
               .drift_content(DriftTest(d$means, w.cov, FALSE)))
  no.plots <- TreeDriftTest(tree, d$means, d$covs, plot = FALSE)
  expect_equal(names(no.plots), names(test.list))
  expect_null(no.plots[[1]]$plot)
})

test_that("TreeDriftTest is independent of the order of means and covariances",
{
  local_seed(4)
  library(ape)
  data(bird.orders)
  tree <- bird.orders
  d <- .drift_tree_data(tree)
  covs.shuffled <- d$covs[sample(length(d$covs))]
  means.shuffled <- d$means[sample(length(d$means))]
  w.cov <- PhyloW(tree, d$covs)$'32'
  test.list <- TreeDriftTest(tree, means.shuffled, covs.shuffled)
  expect_equal(length(test.list), 13)
  expect_equal(.drift_content(test.list[["32"]]),
               .drift_content(DriftTest(d$means[tree$tip.label[getDescendants(tree, 32)[getDescendants(tree, 32) <= 23]]],
                                        w.cov, FALSE)))
})

test_that("TreeDriftTest and PlotTreeDriftTest work with node labels", {
  local_seed(5)
  tree <- ape::read.tree(text = "(((t1:1,t2:1)n1:1,(t3:1,t4:1)n2:1)n3:1,((t5:1,t6:1)n4:1,(t7:1,t8:1)n5:1)n6:1)root;")
  d <- .drift_tree_data(tree)
  res <- TreeDriftTest(tree, d$means, d$covs, plot = FALSE)
  expect_true(all(c("n3", "n6", "root") %in% names(res)))
  expect_false(anyNA(names(res)))
  root.means <- d$means[tree$tip.label]
  expect_equal(.drift_content(res$root),
               .drift_content(DriftTest(root.means, PhyloW(tree, d$covs)$root, FALSE)))
  pdf(NULL); on.exit(dev.off())
  expect_no_error(PlotTreeDriftTest(res, tree))
  # unlabelled trees use node numbers
  tree2 <- tree; tree2$node.label <- NULL
  res2 <- TreeDriftTest(tree2, d$means, d$covs, plot = FALSE)
  expect_true(all(grepl("^[0-9]+$", names(res2))))
  expect_no_error(PlotTreeDriftTest(res2, tree2))
  # repeated labels fall back to numbers with a warning
  tree3 <- tree; tree3$node.label[] <- "100"
  expect_warning(res3 <- TreeDriftTest(tree3, d$means, d$covs, plot = FALSE), "node numbers")
  expect_equal(names(res3), names(res2))
  expect_error(PlotTreeDriftTest(res, tree2), "node names")
})

test_that("TreeDriftTest needs a node with at least 4 tips", {
  tree <- ape::read.tree(text = "((t1:1,t2:1):1,t3:1);")
  d <- .drift_tree_data(tree)
  expect_error(TreeDriftTest(tree, d$means, d$covs), "at least 4")
})

test_that("DriftTest sampling correction works",
{
  local_seed(42)
  p <- 6; m <- 21
  cov.matrix <- RandomMatrix(p, 1, 1, 10, LKJ = FALSE)
  means <- array(rnorm(m * p, sd = 3), c(m, p))
  n <- rep(30, m)
  uncorrected <- DriftTest(means, cov.matrix, FALSE)
  expect_equal(.drift_content(uncorrected), 
               .drift_content(DriftTest(means, cov.matrix, FALSE, sample.sizes = NULL)))
  corrected <- DriftTest(means, cov.matrix, FALSE, sample.sizes = n)
  expect_equal(names(corrected), names(uncorrected))
  expected <- exp(uncorrected$log.between_group_variance) - eigen(cov.matrix)$values / 30
  expected[expected <= 0] <- NA
  expect_equal(corrected$log.between_group_variance, log(expected))
  expect_error(DriftTest(means, cov.matrix, FALSE, sample.sizes = n[-1]))
  # pure sampling error: most axes have no divergence beyond sampling and are dropped
  zero.means <- t(replicate(m, colMeans(MASS::mvrnorm(30, rep(0, p), cov.matrix))))
  expect_warning(try(DriftTest(zero.means, cov.matrix, FALSE, sample.sizes = n), silent = TRUE))
})
