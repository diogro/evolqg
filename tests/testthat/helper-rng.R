# Set a seed, and optionally an old RNG version for tests whose reference
# values were computed under R < 3.6, for the rest of the calling test. The
# previous RNG state comes back when the test ends. Use this instead of
# set.seed() or RNGversion(), which leak into later test files.
local_seed <- function(seed, rng.version = NULL, .local_envir = parent.frame()) {
  if (is.null(rng.version)) {
    withr::local_seed(seed, .local_envir = .local_envir,
                      .rng_kind = "Mersenne-Twister",
                      .rng_normal_kind = "Inversion",
                      .rng_sample_kind = "Rejection")
  } else {
    withr::local_rng_version(rng.version, .local_envir = .local_envir)
    withr::local_seed(seed, .local_envir = .local_envir,
                      .rng_kind = "Mersenne-Twister",
                      .rng_normal_kind = "Inversion",
                      .rng_sample_kind = "Rounding")
  }
}
