.pm_tree <- function() ape::read.tree(text = "(((a:1,b:1):1,(c:1,d:1):1):1,(e:2,f:2):1);")
.pm_matrix <- function(tips, seed) {
  local_seed(seed, .local_envir = parent.frame())
  x <- matrix(rnorm(length(tips)^2), length(tips))
  m <- (x + t(x))/2
  dimnames(m) <- list(tips, tips)
  m
}

test_that("PhyloMantel returns rho and a valid probability", {
  tree <- .pm_tree()
  m1 <- .pm_matrix(tree$tip.label, 1)
  m2 <- .pm_matrix(tree$tip.label, 2)
  local_seed(10)
  out <- PhyloMantel(tree, m1, m2, permutations = 99)
  expect_named(out, c("rho", "Probability"))
  expect_equal(out[["rho"]], cor(m1[lower.tri(m1)], m2[lower.tri(m2)]))
  expect_true(out[["Probability"]] >= 1/100 && out[["Probability"]] <= 1)
  # (count + 1)/(permutations + 1): never exactly zero
  expect_equal(out[["Probability"]] * 100, round(out[["Probability"]] * 100))
})

test_that("PhyloMantel matches matrices to the tips by name", {
  tree <- .pm_tree()
  tips <- tree$tip.label
  m1 <- .pm_matrix(tips, 1)
  m2 <- .pm_matrix(tips, 2)
  ord <- c("f", "c", "a", "e", "b", "d")
  local_seed(10)
  ref <- PhyloMantel(tree, m1, m2, permutations = 30)
  local_seed(10)
  shuffled <- PhyloMantel(tree, m1[ord, ord], m2[rev(ord), rev(ord)], permutations = 30)
  expect_equal(shuffled, ref)
  bad <- m2; colnames(bad)[1] <- "zzz"
  expect_error(PhyloMantel(tree, m1, bad, permutations = 5), "All tip labels")
})

test_that("PhyloMantel uses ComparisonFunc for the null", {
  tree <- .pm_tree()
  m1 <- .pm_matrix(tree$tip.label, 1)
  m2 <- .pm_matrix(tree$tip.label, 2)
  calls <- 0
  const <- function(x, y) { calls <<- calls + 1; 5 }
  local_seed(1)
  out <- PhyloMantel(tree, m1, m2, permutations = 20, ComparisonFunc = const)
  expect_equal(calls, 21)
  expect_equal(out, c(rho = 5, Probability = 1))
  # a function that is always smaller under permutation gives the smallest probability
  perm.small <- local({ first <- TRUE; function(x, y) { if (first) { first <<- FALSE; 10 } else 0 } })
  out <- PhyloMantel(tree, m1, m2, permutations = 20, ComparisonFunc = perm.small)
  expect_equal(out[["Probability"]], 1/21)
})

test_that("PhyloMantel accepts triangular input", {
  tree <- .pm_tree()
  m1 <- .pm_matrix(tree$tip.label, 1)
  m2 <- .pm_matrix(tree$tip.label, 2)
  tri <- m1; tri[upper.tri(tri)] <- NA
  local_seed(3)
  full <- PhyloMantel(tree, m1, m2, permutations = 25)
  local_seed(3)
  from.tri <- PhyloMantel(tree, tri, m2, permutations = 25)
  expect_equal(from.tri, full)
  tri0 <- m1; tri0[upper.tri(tri0)] <- 0
  local_seed(3)
  expect_equal(PhyloMantel(tree, tri0, m2, permutations = 25), full)
  expect_error(PhyloMantel(tree, m1 + upper.tri(m1), m2, permutations = 5), "symmetric or triangular")
})

test_that("PhyloMantel works on the dentus example", {
  data(dentus); data(dentus.tree)
  cor.matrices <- dlply(dentus, .(species), function(x) cor(x[-5]))
  comparisons <- MatrixCor(cor.matrices)
  sp.means <- dlply(dentus, .(species), function(x) colMeans(x[-5]))
  mh.dist <- MultiMahalanobis(means = sp.means, cov.matrix = PhyloW(dentus.tree, cor.matrices)$'6')
  local_seed(1)
  out <- PhyloMantel(dentus.tree, comparisons, mh.dist, k = 10000, permutations = 49)
  expect_named(out, c("rho", "Probability"))
  expect_true(is.finite(out[["rho"]]))
})
