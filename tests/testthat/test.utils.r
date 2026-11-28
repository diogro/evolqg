# Internal replacements for functions from mvtnorm, MCMCpack, coda, vegan and
# expm. The reference values below were computed with mvtnorm 1.4.2,
# MCMCpack, coda, vegan 2.7.3 / permute 0.9.10 and expm, using the inputs
# defined here. The `seed` entries are the next runif() after each call, and
# check that the random number stream advanced exactly as in the original.

S <- matrix(c(2, 0.5, 0.3, 0.5, 1, 0.2, 0.3, 0.2, 1.5), 3)
A <- cov2cor(crossprod(matrix(sin(1:60), 10)))
B <- cov2cor(crossprod(matrix(cos(1:60)^2, 10)))
A9 <- cov2cor(crossprod(matrix(sin(1:180), 20)))
B9 <- cov2cor(crossprod(matrix(cos(1:180)^2, 20)))
A7 <- A9[1:7, 1:7]
B7 <- B9[1:7, 1:7]
ref <- list(rmvnorm_eigen = c(2.8344305634104234, 1.9549687179523989, 
3.3073975352476475, 1.4309631359429114, 1.7630647970249953, 2.5167254896572708, 
2.384500975450381, 3.4443182837643018, 3.5449377985833728, 2.9721788667772464, 
5.6131429250237854, 5.8822769010906688), seed1 = 0.082437558099627495, 
    rmvnorm_svd = c(1.8344305634104232, 0.95496871795239846, 
    2.307397535247647, 0.43096313594291003, -0.23693520297500459, 
    0.51672548965727017, 0.38450097545038148, 1.444318283764302, 
    0.54493779858337232, -0.027821133222753262, 2.6131429250237845, 
    2.8822769010906684), seed2 = 0.082437558099627495, dmvnorm = c(0.040036970013387095, 
    0.0030503324775030031), riwish = c(0.22892074227350351, 0.14035272682449856, 
    -0.022495364749020728, 0.14035272682449856, 0.16568257658144428, 
    -0.012340463088855324, -0.022495364749020728, -0.012340463088855325, 
    0.20890449125455746), seed3 = 0.25542882434092462, hpd = c(-2.6564554209047762, 
    -2.4142076499466323, -1.3682810444192945, 2.0184237138770413, 
    1.214674699172599, 1.5757275197919773), hpd1 = c(-2.6564554209047762, 
    2.0184237138770413), mantel6 = c(-0.4757801498600473, 0.93000000000000005, 
    99), seed4 = 0.77582504274323583, mantel9 = c(-0.26798006792496609, 
    0.96999999999999997, 199), seed5 = 0.22707688226364553, mantel7 = c(-0.38321692382657735, 
    0.97299999999999998, 999), seed6 = 0.81061906577087939, sqrtm = c(1.3947299170887308, 
    0.20729611172432355, 0.10842868827777428, 0.20729611172432352, 
    0.97492073797196543, 0.080980718298713408, 0.10842868827777419, 
    0.080980718298713256, 1.2172449806107997), logm = c(0.63153618376453136, 
    0.35240822525393029, 0.15639575152736324, 0.35240822525393012, 
    -0.091287009050569523, 0.13427507153031001, 0.15639575152736332, 
    0.13427507153030988, 0.38202362883769514))

test_that(".rmvnorm reproduces mvtnorm::rmvnorm and its random stream", {
  local_seed(42)
  x <- .rmvnorm(4, mean = c(a = 1, b = 2, c = 3), sigma = S)
  expect_equal(as.vector(x), ref$rmvnorm_eigen, tolerance = 1e-14)
  expect_equal(colnames(x), c("a", "b", "c"))
  expect_identical(runif(1), ref$seed1)
  local_seed(42)
  x <- .rmvnorm(4, sigma = S, method = "svd")
  expect_equal(as.vector(x), ref$rmvnorm_svd, tolerance = 1e-14)
  expect_identical(runif(1), ref$seed2)
  local_seed(42)
  y <- .rmvnorm(4, sigma = S, method = "svd", root = .mvn_root(S, "svd"))
  expect_identical(x, y)
})

test_that(".dmvnorm reproduces mvtnorm::dmvnorm", {
  x <- rbind(c(0, 0, 0), c(1, -1, 2))
  expect_equal(.dmvnorm(x, sigma = S), ref$dmvnorm, tolerance = 1e-14)
})

test_that(".riwish reproduces MCMCpack::riwish and its random stream", {
  local_seed(42)
  expect_equal(as.vector(.riwish(10, S)), ref$riwish, tolerance = 1e-14)
  expect_identical(runif(1), ref$seed3)
})

test_that(".hpd_interval reproduces coda::HPDinterval", {
  local_seed(42)
  z <- matrix(rnorm(60), 20, dimnames = list(NULL, c("x", "y", "z")))
  h <- .hpd_interval(z, prob = 0.9)
  expect_equal(as.vector(h), ref$hpd)
  expect_equal(dimnames(h), list(c("x", "y", "z"), c("lower", "upper")))
  expect_equal(attr(h, "Probability"), 0.9)
  h1 <- .hpd_interval(z[, 1], prob = 0.9)
  expect_equal(as.vector(h1), ref$hpd1)
  expect_equal(rownames(h1), "var1")
})

test_that(".mantel_test reproduces vegan::mantel and its random stream", {
  # 6 traits: complete enumeration of the 719 permutations, 99 sampled
  local_seed(42)
  m <- .mantel_test(A, B, permutations = 99)
  expect_equal(c(m$statistic, m$signif, m$permutations), ref$mantel6, tolerance = 1e-14)
  expect_identical(runif(1), ref$seed4)
  # 9 traits: permutations drawn one at a time
  local_seed(42)
  m <- .mantel_test(A9, B9, permutations = 199)
  expect_equal(c(m$statistic, m$signif, m$permutations), ref$mantel9, tolerance = 1e-14)
  expect_identical(runif(1), ref$seed5)
  # 7 traits: the 5039 permutations are enumerated, then 999 rows sampled
  local_seed(42)
  m <- .mantel_test(A7, B7, permutations = 999)
  expect_equal(c(m$statistic, m$signif, m$permutations), ref$mantel7, tolerance = 1e-14)
  expect_identical(runif(1), ref$seed6)
})

test_that(".all_permutations lists every permutation once, identity first", {
  P <- .all_permutations(4)
  expect_equal(dim(P), c(24, 4))
  expect_equal(P[1, ], 1:4)
  expect_equal(anyDuplicated(apply(P, 1, paste, collapse = "")), 0)
})

test_that("symmetric matrix functions match expm", {
  expect_equal(as.vector(.sym_sqrtm(S)), ref$sqrtm, tolerance = 1e-12)
  expect_equal(as.vector(.sym_logm(S)), ref$logm, tolerance = 1e-12)
  expect_equal(.sym_expm(.sym_logm(S)), S, tolerance = 1e-12)
  expect_equal(.sym_inv_sqrtm(S) %*% S %*% .sym_inv_sqrtm(S), diag(3), tolerance = 1e-12)
})

test_that(".evq_lapply matches lapply when not parallel", {
  x <- list(a = 1:3, b = 4:6)
  expect_identical(.evq_lapply(x, sum), lapply(x, sum))
  expect_identical(.evq_lapply(x, function(v, k) v * k, k = 2), lapply(x, function(v, k) v * k, k = 2))
})

test_that(".is_cor_matrix tolerates rounding", {
  expect_true(.is_cor_matrix(diag(3) + 1e-12))
  expect_false(.is_cor_matrix(diag(c(1, 2, 1))))
})
