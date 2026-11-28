# the original (quadratic, 32-bit) algorithm, kept as a reference
.old_combine_hypot <- function(modularity.hypot) {
  DecToBin <- function(x) sapply(strsplit(paste(rev(intToBits(x))), ""), `[[`, 2)
  n.hypots <- dim(modularity.hypot)[2]
  if (is.null(colnames(modularity.hypot))) colnames(modularity.hypot) <- 1:n.hypots
  counter <- 2^n.hypots - 1
  hypot_list <- list(null = diag(dim(modularity.hypot)[1]))
  k <- 2
  for (i in seq(counter)) {
    mask <- DecToBin(i)
    mask <- as.logical(as.numeric((mask[(32 - (n.hypots - 1)):32])))
    if (sum(mask) > 1) new_hypot <- CreateHypotMatrix(modularity.hypot[, mask])[[sum(mask) + 1]]
    else new_hypot <- CreateHypotMatrix(modularity.hypot[, mask])
    diag(new_hypot) <- 1
    if (!any(sapply(hypot_list, function(x) all(x == new_hypot)))) {
      hypot_list[[k]] <- new_hypot
      names(hypot_list)[[k]] <- paste(colnames(modularity.hypot)[mask], collapse = "_")
      k <- k + 1
    }
  }
  hypot_list
}

test_that("CombineHypot gives the same matrices and names as the original", {
  local_seed(1)
  hyp <- matrix(sample(c(0, 1), 8*5, TRUE), 8, 5)
  colnames(hyp) <- letters[1:5]
  new <- CombineHypot(hyp)
  masks <- attr(new, "masks")
  attr(new, "masks") <- NULL
  expect_identical(new, .old_combine_hypot(hyp))
  expect_length(masks, length(new))
  expect_false(any(masks[[1]]))
  # the stored masks recover the combination names
  expect_equal(vapply(masks[-1], function(m) paste(colnames(hyp)[m], collapse = "_"), ""),
               names(new)[-1])
  # unnamed columns, with duplicated hypotheses
  hyp2 <- cbind(hyp[, 1:2], hyp[, 1])
  new2 <- CombineHypot(hyp2)
  attr(new2, "masks") <- NULL
  expect_identical(new2, .old_combine_hypot(hyp2))
  # single hypothesis
  expect_equal(CombineHypot(c(1, 1, 0)), `diag<-`(outer(c(1, 1, 0), c(1, 1, 0)), 1))
})

test_that("CombineHypot limits the number of hypotheses", {
  expect_error(CombineHypot(matrix(1, 3, 21)), "Too many hypotheses")
  expect_warning(CombineHypot(matrix(c(1, 0, 0, 1), 2, 2)[, c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2)]),
                 "infeasible")
})

test_that("CreateHypotMatrix returns the expected structure", {
  hyp <- cbind(c(1, 1, 0, 0), c(0, 0, 1, 1))
  out <- CreateHypotMatrix(hyp)
  expect_length(out, 3)
  expect_equal(out[[3]], matrix(c(1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1), 4),
               ignore_attr = TRUE)
  expect_equal(CreateHypotMatrix(c(1, 0, 1)), outer(c(1, 0, 1), c(1, 0, 1)))
})

test_that("Partition2HypotMatrix always returns a matrix", {
  expect_equal(Partition2HypotMatrix(c(1, 2, 1, 3)),
               cbind(c(1, 0, 1, 0), c(0, 1, 0, 0), c(0, 0, 0, 1)))
  one <- Partition2HypotMatrix(c(1, 1, 1))
  expect_true(is.matrix(one))
  expect_equal(dim(one), c(3, 1))
  single <- Partition2HypotMatrix(5)
  expect_equal(dim(single), c(1, 1))
  expect_equal(colnames(Partition2HypotMatrix(c("a", "b", "a"))), c("a", "b"))
  expect_error(Partition2HypotMatrix(c(1, NA, 2)), "NA")
})
