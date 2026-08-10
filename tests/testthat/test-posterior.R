## Statistical checks. These run long chains and are skipped on CRAN.

data(dem2gbp, package = "bayesGARCH")

pmean <- function(M, l.bi) {
  X <- do.call(rbind, lapply(M, function(z) as.matrix(z)[-seq_len(l.bi), , drop = FALSE]))
  colMeans(X)
}

test_that("the Normal special case reproduces the dem2gbp benchmark", {
  skip_on_cran()
  ## lambda = 100, delta = 500 centres nu around 500, i.e. approximately Normal
  ## innovations. The GARCH(1,1)-N fit of dem2gbp is an established benchmark:
  ## alpha0 = 0.0108, alpha1 = 0.153, beta = 0.806 (Bollerslev and Ghysels, 1996;
  ## McCullough and Renfro, 1999; Brooks, Burke and Persand, 2001).
  set.seed(7)
  M <- bayesGARCH(dem2gbp, lambda = 100, delta = 500,
                  control = list(n.chain = 2, l.chain = 10000, refresh = 0))
  m <- pmean(M, 2000)

  expect_equal(unname(m["alpha0"]), 0.0108, tolerance = 0.15)
  expect_equal(unname(m["alpha1"]), 0.153,  tolerance = 0.15)
  expect_equal(unname(m["beta"]),   0.806,  tolerance = 0.10)
})

test_that("addPriorConditions targets the constrained posterior", {
  skip_on_cran()
  ## Reference values for p(psi | y) restricted to A = {alpha1 + beta < 1} on the
  ## full dem2gbp series, obtained from an independent adaptive random-walk
  ## Metropolis sampler sharing no code with this package (four over-dispersed
  ## chains, 150k draws each, Rhat = 1.0000, ESS > 12000):
  ##
  ##   alpha0 = 0.00513   alpha1 = 0.1427   beta = 0.8484   nu = 4.566
  ##
  ## They are reproduced independently by conditioning: the unconstrained draws
  ## that happen to satisfy alpha1 + beta < 1 are an exact sample from the same
  ## distribution, and give 0.00515 / 0.1426 / 0.8485 / 4.543.
  ##
  ## Up to version 2.1.10 the constraint was applied to the whole block and the
  ## sampler converged instead to 0.00841 / 0.1930 / 0.7953 / 4.494 -- alpha0
  ## overstated by roughly 60% and alpha1 by roughly 35%. The tolerances below
  ## are wide enough for Monte Carlo error at this chain length and far tighter
  ## than that bias.
  ref <- c(alpha0 = 0.00513, alpha1 = 0.1427, beta = 0.8484, nu = 4.566)

  sv <- matrix(c(0.002, 0.08, 0.90, 3,
                 0.012, 0.25, 0.70, 6), 2, 4, byrow = TRUE)
  set.seed(5)
  M <- bayesGARCH(dem2gbp, control = list(
    l.chain = 30000, start.val = sv, refresh = 0,
    addPriorConditions = function(psi) psi[2] + psi[3] < 1))
  m <- pmean(M, 5000)

  ## the constraint itself is never violated
  expect_true(all(vapply(M, function(z) all(z[, 2] + z[, 3] < 1), NA)))

  expect_equal(unname(m["alpha0"]), unname(ref["alpha0"]), tolerance = 0.15)
  expect_equal(unname(m["alpha1"]), unname(ref["alpha1"]), tolerance = 0.10)
  expect_equal(unname(m["beta"]),   unname(ref["beta"]),   tolerance = 0.05)
  expect_equal(unname(m["nu"]),     unname(ref["nu"]),     tolerance = 0.10)
})

test_that("an inactive constraint leaves the sampler unchanged", {
  skip_on_cran()
  ## addPriorConditions that is always TRUE must reproduce the unconstrained run
  ## draw for draw.
  y <- dem2gbp[1:400]
  set.seed(31)
  a <- bayesGARCH(y, control = list(l.chain = 300, refresh = 0))
  set.seed(31)
  b <- bayesGARCH(y, control = list(l.chain = 300, refresh = 0,
                                    addPriorConditions = function(psi) TRUE))
  expect_equal(as.matrix(a[[1]]), as.matrix(b[[1]]))
})
