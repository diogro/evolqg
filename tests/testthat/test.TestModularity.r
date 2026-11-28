test_that("TestModularity returns correct results",
          {
            local_seed(43, "3.5.0")
            cor.matrix <- RandomMatrix(10)
            rand.hipots <- matrix(sample(c(1, 0), 30, replace=TRUE), 10, 3)
            hip.array <- CreateHypotMatrix(rand.hipots)
            expect_equal(sum(laply(hip.array, isSymmetric)), 4)
            expect_type(hip.array, "list")
            expect_equal(dim(hip.array[[1]]), c(10, 10))
            mod.test <- TestModularity(cor.matrix, rand.hipots)
            expect_equal(dim(mod.test), c(4, 6))
            expect_equal(mod.test[,"hypothesis"], c("1", "2", "3", "Full Integration"))
            expect_equal(colnames(mod.test), c("hypothesis", "rho", "Probability", "AVG+", "AVG-", "AVG Ratio"))
            expect_true(all(mod.test[1:4, 2:5] >= -1))
            expect_true(all(mod.test[1:4, 2:5] <= 1))
            expect_equal(mod.test[,4]/mod.test[,5], mod.test[,6])
          }
)

test_that("TestModularity returns correct results for Modularity Hypothesis Index",
          {
            local_seed(43, "3.5.0")
            cov.matrix <- RemoveSize(RandomMatrix(11))
            rand.hipots <- matrix(sample(c(1, 0), 33, replace=TRUE), 11, 3)
            hip.array <- CreateHypotMatrix(rand.hipots)
            mod.test <- TestModularity(cov.matrix, rand.hipots, MHI = TRUE)
            expect_equal(dim(mod.test), c(4, 6))
            expect_equal(mod.test[,"hypothesis"], c("1", "2", "3", "Full Integration"))
            expect_equal(colnames(mod.test), c("hypothesis", "rho", 
                                                "Probability", "AVG+", "AVG-", "MHI"))
            expect_true(all(mod.test[1:4, 2:5] >= -1))
            expect_true(all(mod.test[1:4, 2:5] <= 1))
            expect_equal((mod.test[,4] - mod.test[,5])/CalcEigenVar(cov.matrix,sd = TRUE, rel = FALSE), mod.test[,6])
          }
)

test_that("MantelModTest returns correct results for Modularity Hypothesis Index",
          {
            # Create a single modularity hypothesis:
            hypot = rep(c(1, 0), each = 6)
            cor.hypot = CreateHypotMatrix(hypot)
            
            # First with an unstructured matrix:
            expect = c(rho = -0.0393920752173474, Probability = 0.556443556443556, 
`AVG+` = -0.0596580385711818, `AVG-` = -0.0271784073214975, MHI = -0.0254288119834227)
            local_seed(42, "3.5.0")
            un.cor = RandomMatrix(12, LKJ = FALSE)
            result = MantelModTest(cor.hypot, RemoveSize(un.cor), MHI = TRUE)
            expect_equal(result, expect)
            
            # Now with a modular matrix:
            expect = c(rho = 1, Probability = 0.001998001998002, `AVG+` = 0.8, 
                       `AVG-` = 0.3, MHI = 0.325128044)
            local_seed(42, "3.5.0")
            hypot.mask = matrix(as.logical(cor.hypot), 12, 12)
            mod.cor = matrix(NA, 12, 12)
            mod.cor[ hypot.mask] = 0.8 # within-modules
            mod.cor[!hypot.mask] = 0.3 # between-modules
            diag(mod.cor) = 1
            result = MantelModTest(cor.hypot, mod.cor, MHI = TRUE)
            expect_equal(result, expect)
            avg = CalcAVG(cor.hypot, mod.cor, MHI = TRUE)
            expect_equal(avg, result[3:5])
          }
)

test_that("MantelModTest returns correct results for non-landmark data",
          {
            # Create a single modularity hypothesis:
            hypot = rep(c(1, 0), each = 6)
            cor.hypot = CreateHypotMatrix(hypot)
             
            # First with an unstructured matrix:
            expect = structure(c(0.0810901974853795, 0.275724275724276, 0.158055575956429, 0.061087664745461, 2.587356655636), names = c("rho", "Probability", "AVG+", "AVG-", "AVG Ratio")) 
            local_seed(42, "3.5.0")
            un.cor = RandomMatrix(12, LKJ = FALSE)
            result = MantelModTest(cor.hypot, un.cor)
            expect_equal(result, expect)
            
            # Now with a modular matrix:
            expect = structure(c(0.99102012957595, 0.001998001998002, 0.87040471785857, 0.354135870234584, 2.45782704045767), names = c("rho", "Probability", "AVG+", "AVG-", "AVG Ratio"))
            local_seed(42, "3.5.0")
            hypot.mask = matrix(as.logical(cor.hypot), 12, 12)
            mod.cor = matrix(NA, 12, 12)
            mod.cor[ hypot.mask] = runif(length(mod.cor[ hypot.mask]), 0.8, 0.9) # within-modules
            mod.cor[!hypot.mask] = runif(length(mod.cor[!hypot.mask]), 0.3, 0.4) # between-modules
            # Replicate the lower triangle to make it symmetric:
            mod.cor[upper.tri(mod.cor)] = t(mod.cor)[upper.tri(mod.cor)]
            diag(mod.cor) = 1
            
            result = MantelModTest(cor.hypot, mod.cor)
            expect_equal(result, expect)
            result = MantelModTest(cor.hypot, mod.cor, MHI = TRUE)
            expect_equal(CalcAVG(cor.hypot, mod.cor), result[3:5])
          })

test_that("MantelModTest returns correct results for landmark data",
          {
            # Create a single modularity hypothesis:
            hypot = rep(c(1, 0), each = 6)
            cor.hypot = CreateHypotMatrix(hypot)
            
            # First with an unstructured matrix:
            expect = structure(c(0.112390456689369, 310/1001, 0.209791381269483, 0.0703378394237725, 2.98262475771438, -0.0679013877016532), names = c("rho", "Probability", "AVG+", "AVG-", "AVG Ratio", "AVG within landmark"))
            local_seed(42, "3.5.0")
            un.cor = RandomMatrix(12, LKJ = FALSE)
            result = MantelModTest(cor.hypot, un.cor, landmark.dim = 2)
            expect_equal(result, expect)
            
            # Now with a modular matrix:
            expect = structure(c(1, 0, 0.8, 0.3, 2.66666666666667, 0.1), names = c("rho","Probability", "AVG+", "AVG-", "AVG Ratio", "AVG within landmark"))
            hypot.mask = matrix(as.logical(cor.hypot), 12, 12)
            mod.cor = matrix(NA, 12, 12)
            mod.cor[ hypot.mask] = 0.8 # within-modules
            mod.cor[!hypot.mask] = 0.3 # between-modules
            mod.cor[evolqg:::CreateWithinLandMat(6, 2)] = 0.1# within-land
            diag(mod.cor) = 1
            result = MantelModTest(cor.hypot, mod.cor, landmark.dim = 2)
            # rho = 1 is reached by every landmark permutation that keeps the
            # module together (3! * 3! of 6! = 5%), so ties make the p-value ~0.05
            expect_equal(result[-2], expect[-2])
            expect_equal(unname(result["Probability"]), 0.05, tolerance = 0.3)
            expect_equal(CalcAVG(cor.hypot, mod.cor, landmark.dim = 2, MHI = FALSE), result[3:6])
            
            mod.cor = matrix(NA, 12, 12)
            mod.cor[ hypot.mask] = 0.8 # within-modules
            mod.cor[!hypot.mask] = 0.3 # between-modules
            mod.cor[evolqg:::CreateWithinLandMat(4, 3)] = 0.1 # within-land
            diag(mod.cor) = 1
            result = MantelModTest(cor.hypot, mod.cor, landmark.dim = 3)
            # here 2! * 2! of 4! landmark permutations (1/6) keep the module together
            expect_equal(result[-2], expect[-2])
            expect_equal(unname(result["Probability"]), 1/6, tolerance = 0.2)
            expect_equal(CalcAVG(cor.hypot, mod.cor, landmark.dim = 3, MHI = FALSE), result[3:6])
            })

test_that("MantelModTest trows errors",
          {
            hypot = rep(c(1, 0), each = 6)
            cor.hypot = CreateHypotMatrix(hypot)
            un.cor = RandomMatrix(12)
            expect_error(MantelModTest(un.cor, cor.hypot), 
                         "modularity hypothesis matrix should be binary")
            expect_error(MantelModTest(cor.hypot, cor.hypot, landmark.dim = 4), 
                         "landmark.dim should be either 2 or 3 dimensions")
            expect_error(CalcAVG(un.cor, cor.hypot), 
                         "modularity hypothesis matrix should be binary")
            expect_error(CalcAVG(cor.hypot, cor.hypot, landmark.dim = 4), 
                         "landmark.dim should be either 2 or 3 dimensions")
          })

test_that("TestModularity shares permutations and keeps the output structure", {
  local_seed(7)
  cor.matrix <- RandomMatrix(10)
  hyp <- matrix(c(rep(c(1, 0), each = 5), rep(c(0, 1), each = 5), rep(c(1, 0), 5)), 10, 3)
  colnames(hyp) <- c("A", "B", "C")
  local_seed(1)
  out <- TestModularity(cor.matrix, hyp, permutations = 100)
  expect_s3_class(out, "data.frame")
  expect_equal(names(out), c("hypothesis", "rho", "Probability", "AVG+", "AVG-", "AVG Ratio"))
  expect_equal(out$hypothesis, c("A", "B", "C", "Full Integration"))
  # rho, AVG+ and AVG- do not depend on the permutations
  single <- lapply(CreateHypotMatrix(hyp), function(h) MantelModTest(h, cor.matrix, permutations = 100))
  expect_equal(out$rho, unname(vapply(single, `[[`, 0, "rho")))
  expect_equal(out[["AVG Ratio"]], unname(vapply(single, `[[`, 0, "AVG Ratio")))
  expect_true(all(out$Probability >= 1/101 & out$Probability <= 1))
  # hypotheses given the same permutations get the same p-value if they give the same statistic
  same <- TestModularity(cor.matrix, cbind(hyp[, 1], hyp[, 1]), permutations = 100)
  expect_equal(same$Probability[1], same$Probability[2])
})

test_that("MantelModTest passes withinLandmark to MantelCor", {
  local_seed(3)
  cor.hypot <- CreateHypotMatrix(rep(c(1, 0), each = 6))
  un.cor <- RandomMatrix(12, LKJ = FALSE)
  local_seed(1)
  without <- MantelModTest(cor.hypot, un.cor, landmark.dim = 2, permutations = 20)
  local_seed(1)
  within <- MantelModTest(cor.hypot, un.cor, landmark.dim = 2, withinLandmark = TRUE, permutations = 20)
  expect_false(isTRUE(all.equal(without[["rho"]], within[["rho"]])))
  # also through the list method and TestModularity
  hyp <- cbind(rep(c(1, 0), each = 6), rep(c(0, 1), each = 6))
  a <- TestModularity(un.cor, hyp, permutations = 20, landmark.dim = 2)
  b <- TestModularity(un.cor, hyp, permutations = 20, landmark.dim = 2, withinLandmark = TRUE)
  expect_false(isTRUE(all.equal(a$rho, b$rho)))
})

test_that("MantelModTest checks the hypothesis", {
  hypot <- CreateHypotMatrix(rep(c(1, 0), each = 6))
  cm <- RandomMatrix(12)
  bad <- hypot; bad[1, 2] <- bad[2, 1] <- 2 # still has 0s and 1s
  expect_error(MantelModTest(bad, cm), "should be binary")
  expect_warning(MantelModTest(matrix(0, 12, 12), cm, permutations = 10), "no within-module")
  expect_warning(MantelModTest(`diag<-`(matrix(0, 12, 12), 1), cm, permutations = 10), "no within-module")
  expect_warning(MantelModTest(matrix(1, 12, 12), cm, permutations = 10), "single module")
  expect_warning(MantelModTest(hypot, cm, permutations = 10), NA)
})

test_that("MantelModTest.list forwards arguments", {
  local_seed(4)
  cm <- RandomMatrix(12)
  hl <- CreateHypotMatrix(cbind(rep(c(1, 0), each = 6), rep(c(0, 1), each = 6)))
  out <- MantelModTest(hl, cm, permutations = 20, MHI = TRUE)
  expect_true("MHI" %in% names(out))
  expect_true(".id" %in% names(out))
  # extra arguments reach MantelCor
  expect_error(MantelModTest(hl, cm, permutations = 20, landmark.dim = 4), "landmark.dim")
})
