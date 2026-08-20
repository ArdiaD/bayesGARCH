## Fast checks that the sampler stays inside the region it is meant to explore.
## These deliberately do not call skip_on_cran(): up to 2.2.0 the only tests of
## the constrained sampler were the long ones, so a chain that never entered the
## constrained support passed the CRAN suite unnoticed.

data(dem2gbp, package = "bayesGARCH")
y <- dem2gbp[1:250]
csc <- function(psi) psi[2] + psi[3] < 1

test_that("inadmissible starting values are rejected, not silently returned", {
  ## outside the user's own constraint: alpha1 + beta = 1.3
  expect_error(
    bayesGARCH(y, control = list(l.chain = 5, refresh = 0,
                                 start.val = c(0.01, 0.6, 0.7, 5),
                                 addPriorConditions = csc)),
    "does not satisfy 'control\\$addPriorConditions'")

  ## outside the built-in prior support
  expect_error(bayesGARCH(y, control = list(l.chain = 5, start.val = c(0, 0.1, 0.7, 5))),
               "must be positive")
  expect_error(bayesGARCH(y, control = list(l.chain = 5, start.val = c(0.01, -0.1, 0.7, 5))),
               "must be positive")
  expect_error(bayesGARCH(y, control = list(l.chain = 5, start.val = c(0.01, 0.1, 0.7, 2))),
               "'nu' must be greater than 'delta'")
  expect_error(bayesGARCH(y, lambda = 100, delta = 500,
                          control = list(l.chain = 5, start.val = c(0.01, 0.1, 0.7, 100))),
               "'nu' must be greater than 'delta'")
  expect_error(bayesGARCH(y, control = list(l.chain = 5, start.val = c(0.01, 0.1, NA, 5))),
               "non-finite")

  ## the offending row is named
  sv <- matrix(c(0.01, 0.10, 0.70, 5,
                 0.01, 0.60, 0.70, 5), 2, 4, byrow = TRUE)
  expect_error(bayesGARCH(y, control = list(l.chain = 5, start.val = sv,
                                            addPriorConditions = csc)),
               "row 2")
})

test_that("the default starting values are admissible for both innovation settings", {
  ## nu starts at the prior mean delta + 1/lambda, so the Normal special case
  ## (delta = 500) starts inside the prior support rather than below it
  expect_silent(bayesGARCH(y, control = list(l.chain = 3, refresh = 0)))
  expect_silent(bayesGARCH(y, lambda = 100, delta = 500,
                           control = list(l.chain = 3, refresh = 0)))
  m <- bayesGARCH(y, lambda = 100, delta = 500, control = list(l.chain = 3, refresh = 0))
  expect_gt(as.matrix(m[[1]])[1, "nu"], 500)
})

test_that("every returned draw satisfies the requested constraint", {
  set.seed(42)
  M <- bayesGARCH(y, control = list(n.chain = 2, l.chain = 40, refresh = 0,
                                    start.val = c(0.02, 0.10, 0.80, 6),
                                    addPriorConditions = csc))
  for (ch in M) {
    X <- as.matrix(ch)
    expect_true(all(X[, "alpha1"] + X[, "beta"] < 1))
    expect_true(all(X > 0))
    expect_true(all(is.finite(X)))
  }
})

test_that("a constraint no value of nu can satisfy does not exhaust the proposal loop", {
  ## alpha1 + beta < 1 does not involve nu, so no proposal can satisfy it from an
  ## inadmissible state. The constraint was once tested inside the accept-reject
  ## loop, which then ran all 50000 proposals. The draw is now taken from the full
  ## conditional first and the indicator applied once, so the constrained call
  ## consumes exactly the same random numbers as the unconstrained one. Comparing
  ## the state of the generator afterwards tests that deterministically, rather
  ## than by wall-clock time, which is unreliable on a loaded machine.
  fn.nu.full <- bayesGARCH:::fn.nu.full
  n <- 250
  a <- list(y[1:n], c(0.01, 0.6), 0.7, 5, rep(1, n), c(0, 0),
            solve(1000 * diag(1, 2)), 0, 1/1000, 0.01, 2)

  set.seed(11); v1 <- do.call(fn.nu.full, c(a, list(function(psi) TRUE)))
  r1 <- runif(1)
  set.seed(11); v2 <- do.call(fn.nu.full, c(a, list(csc)))
  r2 <- runif(1)

  expect_identical(r1, r2)          # same number of variates drawn
  expect_identical(v2, 5)           # the inadmissible draw is rejected
  expect_false(identical(v1, v2))   # the unconstrained call keeps it
})

test_that("the sampler is equivariant to a power-of-two rescaling of y", {
  ## y is rescaled internally to unit root mean square, exactly so for a
  ## power-of-two factor. Under y -> y*s only alpha0 carries the scale, so with
  ## the prior on alpha0, the starting value and the constraint all expressed in
  ## the new units, the two runs must agree draw for draw.
  ##
  ## The prior has to be transformed as well: it is stated in absolute units, so
  ## holding it fixed while rescaling y genuinely changes the posterior rather
  ## than revealing an asymmetry in the sampler.
  y1 <- dem2gbp[1:400]
  s  <- 1024
  d  <- c(s^2, 1)
  set.seed(21)
  A <- as.matrix(bayesGARCH(y1, control = list(l.chain = 40, refresh = 0,
                                               start.val = c(0.05, 0.1, 0.7, 10)))[[1]])
  set.seed(21)
  B <- as.matrix(bayesGARCH(y1 * s,
                            mu.alpha    = c(0, 0) * d,
                            Sigma.alpha = 1000 * diag(1, 2) * outer(d, d),
                            control = list(l.chain = 40, refresh = 0,
                                           start.val = c(0.05 * s^2, 0.1, 0.7, 10)))[[1]])
  expect_identical(A[, "alpha1"], B[, "alpha1"])
  expect_identical(A[, "beta"],   B[, "beta"])
  expect_identical(A[, "nu"],     B[, "nu"])
  expect_identical(A[, "alpha0"] * s^2, B[, "alpha0"])
})

test_that("badly scaled data no longer stalls or crashes", {
  ## dem2gbp is quoted in percent; a user supplying decimal log-returns passes
  ## y/100. That used to reject every alpha and beta move -- move rates of
  ## exactly zero -- and then abort in the nu update once nu had been dragged to
  ## delta and the latent variances had collapsed.
  for (s in c(1, 1e-2, 1e-6)) {
    set.seed(5)
    M <- bayesGARCH(y * s, control = list(l.chain = 60, refresh = 0))
    mr <- attr(M[[1]], "move.rates")
    expect_gt(mr[["alpha0"]], 0)
    expect_gt(mr[["alpha1"]], 0)
    expect_gt(mr[["beta"]],   0)
    expect_true(all(is.finite(as.matrix(M[[1]]))))
  }
  ## a scale that cannot be represented is refused with a clear message
  expect_error(bayesGARCH(y * 1e-200, control = list(l.chain = 3)), "extreme scale")
  expect_error(bayesGARCH(y * 1e155,  control = list(l.chain = 3)), "extreme scale")
})

test_that("the nu root bracket widens instead of failing", {
  ## a sharply concentrated prior on nu pushes the rejection-envelope rate above
  ## the former fixed upper limit of 500
  expect_silent(bayesGARCH(y[1:150], lambda = 500, delta = 500,
                           control = list(l.chain = 3, refresh = 0)))
  expect_silent(bayesGARCH(y[1:150], lambda = 5000, delta = 500,
                           control = list(l.chain = 3, refresh = 0)))
})
