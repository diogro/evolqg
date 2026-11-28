.tiny_tree <- function() ape::read.tree(text = "((a:1,b:1)n1:1,(c:1,d:1)n2:1)root;")
.tiny_mats <- function() lapply(c(a = 1, b = 5, c = 10, d = 20), function(x) matrix(x))

test_that("PhyloW gives correct results", {
  TREE <- "iris((setosa:1,versicolor:1):1,virginica:2);"
  tree.iris <- ape::read.tree(text = TREE)
  iris.cov.list <- dlply(iris, "Species", function(x) cov(x[,1:4]))
  cov.matrices <- PhyloW(tree.iris, iris.cov.list)
  w_matrix <- CalculateMatrix(lm(as.matrix(iris[,1:4])~iris[,5]))
  expect_equal(cov.matrices$'4', w_matrix)
  versicolor_setosa <- Reduce("+", cov.matrices[1:2])/2
  expect_equal(cov.matrices$'5', versicolor_setosa)
  expect_error(PhyloW(tree.iris, iris.cov.list[1:2]), "All tip labels must be in stat list.")
})

test_that("PhyloW uses node labels and tip order of the data", {
  tree <- .tiny_tree()
  mats <- .tiny_mats()[c("c", "a", "d", "b")] # not in tip order
  out <- PhyloW(tree, mats)
  expect_equal(names(out), c("c", "a", "d", "b", "n2", "n1", "root")) # internal nodes in postorder
  expect_equal(out$n1[1, 1], 3)
  expect_equal(out$n2[1, 1], 15)
  expect_equal(out$root[1, 1], 9)
})

test_that("PhyloW keeps names of tip.sample.size", {
  tree <- .tiny_tree()
  mats <- .tiny_mats()[c("c", "a", "d", "b")]
  sizes <- c(d = 1, b = 3, a = 1, c = 1) # named, in another order
  out <- PhyloW(tree, mats, sizes)
  expect_equal(out$n1[1, 1], (1*1 + 5*3)/4)
  expect_equal(out$n2[1, 1], 15)
  expect_equal(out$root[1, 1], (4*(1*1 + 5*3)/4 + 2*15)/6)
  # unnamed sizes follow the order of tip.data, as before
  out2 <- PhyloW(tree, mats, c(1, 1, 1, 3))
  expect_equal(out2$n1[1, 1], (1*1 + 5*3)/4)
  expect_error(PhyloW(tree, mats, c(a = 1, b = 2)), "tip.sample.size")
})

test_that("PhyloW handles repeated and empty node labels", {
  mats <- .tiny_mats()
  dup <- ape::read.tree(text = "((a:1,b:1)100:1,(c:1,d:1)100:1)100;")
  expect_warning(out <- PhyloW(dup, mats), "node numbers")
  expect_equal(names(out), c("a", "b", "c", "d", "7", "6", "5"))
  expect_equal(out[["6"]][1, 1], 3)
  expect_equal(out[["7"]][1, 1], 15)
  empty <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1)n2:1);")
  expect_warning(out <- PhyloW(empty, mats), "node numbers")
  expect_equal(out[["6"]][1, 1], 3)
  expect_equal(out[["7"]][1, 1], 15)
  nolab <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  expect_warning(out <- PhyloW(nolab, mats), NA)
  expect_equal(names(out), c("a", "b", "c", "d", "7", "6", "5"))
})
