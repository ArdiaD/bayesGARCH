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

test_that("a constraint no value of nu can satisfy does not spin the nu step", {
  ## alpha1 + beta < 1 does not involve nu. Up to 2.2.0 the constraint was tested
  ## inside the accept-reject loop, so an inadmissible state burned all 50000
  ## proposals per iteration. Starting admissibly, the run must stay quick.
  set.seed(1)
  el <- system.time(bayesGARCH(y, control = list(l.chain = 25, refresh = 0,
                                                 start.val = c(0.02, 0.10, 0.80, 6),
                                                 addPriorConditions = csc)))[["elapsed"]]
  expect_lt(el, 5)
})
