## Argument checking and the control list.

data(dem2gbp, package = "bayesGARCH")
y <- dem2gbp[1:300]

test_that("bayesGARCH rejects malformed arguments", {
  expect_error(bayesGARCH(), "'y' is missing")
  expect_error(bayesGARCH(matrix(1:4, 2)), "'y' must be a numeric vector")
  expect_error(bayesGARCH(letters[1:5]), "'y' must be a numeric vector")
  expect_error(bayesGARCH(1), "longer time series")
  expect_error(bayesGARCH(c(1, NA, 2)), "'NA'")
  expect_error(bayesGARCH(y, mu.alpha = 1:3), "appropriate size")
  expect_error(bayesGARCH(y, Sigma.alpha = diag(1, 3)), "appropriate size")
  expect_error(bayesGARCH(y, Sigma.alpha = matrix(c(1, 2, 3, 1), 2)), "not symmetric")
  expect_error(bayesGARCH(y, Sigma.alpha = matrix(c(1, 2, 2, 1), 2)), "positive definite")
  expect_error(bayesGARCH(y, Sigma.beta = -1), "finite positive value")
  expect_error(bayesGARCH(y, lambda = -1), "finite positive value")
  expect_error(bayesGARCH(y, delta = 1), "lower than '2'")
  expect_error(bayesGARCH(y, control = 1), "'control' must be a list")
  expect_error(bayesGARCH(y, control = list(addPriorConditions = 1)), "NULL or a function")
})

test_that("the control list is validated", {
  expect_error(bayesGARCH(y, control = list(l.chain = 1)), "'control\\$l.chain'")
  expect_error(bayesGARCH(y, control = list(l.chain = 10, n.chain = 0)), "'control\\$n.chain'")
  expect_error(bayesGARCH(y, control = list(l.chain = 10, start.val = c(0.1, 0.1, 0.7))),
               "four columns")
})

test_that("start.val and n.chain combine as documented", {
  sv <- matrix(c(0.01, 0.10, 0.70, 10,
                 0.05, 0.30, 0.50, 30), 2, 4, byrow = TRUE)

  ## several rows, n.chain omitted: one chain per row
  a <- bayesGARCH(y, control = list(l.chain = 5, start.val = sv, refresh = 0))
  expect_length(a, 2L)
  expect_equal(unname(vapply(a, function(z) z[1, 1], 0)), c(0.01, 0.05))

  ## several rows AND n.chain: used as they stand (regression test, #2.2.0)
  b <- bayesGARCH(y, control = list(n.chain = 2, l.chain = 5, start.val = sv, refresh = 0))
  expect_length(b, 2L)
  expect_equal(unname(vapply(b, function(z) z[1, 1], 0)), c(0.01, 0.05))

  ## a single row is replicated across chains
  d <- bayesGARCH(y, control = list(n.chain = 3, l.chain = 5,
                                    start.val = c(0.01, 0.1, 0.7, 10), refresh = 0))
  expect_length(d, 3L)
  expect_equal(unname(vapply(d, function(z) z[1, 1], 0)), rep(0.01, 3))

  ## an inconsistent request is reported clearly
  expect_error(bayesGARCH(y, control = list(n.chain = 3, l.chain = 5, start.val = sv)),
               "does not match the number of rows")
})

test_that("the output has the documented shape and is reproducible", {
  set.seed(77)
  a <- bayesGARCH(y, control = list(n.chain = 2, l.chain = 20, refresh = 0))
  set.seed(77)
  b <- bayesGARCH(y, control = list(n.chain = 2, l.chain = 20, refresh = 0))

  expect_s3_class(a, "mcmc.list")
  expect_s3_class(a[[1]], "mcmc")
  expect_named(a, c("chain1", "chain2"))
  expect_equal(dim(a[[1]]), c(20L, 4L))
  expect_equal(colnames(a[[1]]), c("alpha0", "alpha1", "beta", "nu"))
  expect_true(all(vapply(a, function(z) all(z > 0), NA)))
  expect_identical(a, b)

  ## per-chain move rates are reported alongside the draws
  mr <- attr(a[[1]], "move.rates")
  expect_named(mr, c("alpha0", "alpha1", "beta", "nu"))
  expect_true(all(mr >= 0 & mr <= 1))
})

test_that("a 'ts' series is accepted", {
  set.seed(78)
  a <- bayesGARCH(as.numeric(y), control = list(l.chain = 15, refresh = 0))
  set.seed(78)
  b <- bayesGARCH(stats::ts(y, frequency = 12), control = list(l.chain = 15, refresh = 0))
  expect_equal(as.matrix(a[[1]]), as.matrix(b[[1]]))
})

test_that("progress reporting goes through message() and is suppressible", {
  msgs <- capture_messages(bayesGARCH(y, control = list(l.chain = 4, refresh = 1)))
  expect_true(any(grepl("chain: 1", msgs, fixed = TRUE)))
  expect_silent(bayesGARCH(y, control = list(l.chain = 4, refresh = 0)))
})

test_that("formSmpl validates its arguments and batches correctly", {
  M <- matrix(seq_len(40), 10, 4,
              dimnames = list(NULL, c("alpha0", "alpha1", "beta", "nu")))

  expect_error(formSmpl(), "'MCMC' is missing")
  expect_error(formSmpl("a"), "neither a list or a matrix")
  expect_error(formSmpl(M, l.bi = -1), "non-negative whole number")
  expect_error(formSmpl(M, l.bi = 10), "l.bi >= l.chain")
  expect_error(formSmpl(M, batch.size = 0), "larger or equal than 1")
  ## validated 'l.bi' instead of 'batch.size' before 2.2.0
  expect_error(formSmpl(M, batch.size = c(1, 2)), "'batch.size' must be a scalar")

  expect_message(formSmpl(M, l.bi = 2, batch.size = 2), "smpl size: 4")
  r <- suppressMessages(formSmpl(M, l.bi = 2, batch.size = 2))
  expect_s3_class(r, "mcmc")
  expect_equal(nrow(r), 4L)          # rows 3,5,7,9
  expect_equal(as.numeric(r[, 1]), c(3, 5, 7, 9))
})

test_that("unknown control components are rejected", {
  ## 'hypers' is internal and used to bypass the validated hyper-parameters
  expect_error(
    bayesGARCH(y, control = list(l.chain = 3, refresh = 0,
      hypers = list(mu.alpha = c(0, 0), iv.alpha = solve(1000 * diag(1, 2)),
                    mu.beta = 0, iv.beta = 1/1000, c.nu = -5, d.nu = 2))),
    "unknown component")
  ## a mistyped component would otherwise be ignored and the default used
  expect_error(bayesGARCH(y, control = list(l.chian = 5, refresh = 0)), "l.chian")
  expect_silent(bayesGARCH(y, control = list(l.chain = 3, refresh = 0)))
})

test_that("count-like controls must be finite whole numbers", {
  for (bad in list(2.5, NA, NaN, Inf, -1))
    expect_error(bayesGARCH(y, control = list(l.chain = bad)), "control\\$l.chain")
  for (bad in list(1.5, NA, Inf, 0))
    expect_error(bayesGARCH(y, control = list(l.chain = 3, n.chain = bad)), "control\\$n.chain")
  for (bad in list(NA, Inf, -1, 0.5))
    expect_error(bayesGARCH(y, control = list(l.chain = 3, refresh = bad)), "control\\$refresh")
  for (bad in list(NA, Inf, -1, 0.5))
    expect_error(bayesGARCH(y, control = list(l.chain = 3, digits = bad)), "control\\$digits")
})

test_that("non-finite data and hyper-parameters are rejected up front", {
  expect_error(bayesGARCH(c(y[1:50], Inf), control = list(l.chain = 3)), "non-finite")
  expect_error(bayesGARCH(y, mu.beta = Inf, control = list(l.chain = 3)), "'mu.beta' must be finite")
  expect_error(bayesGARCH(y, Sigma.beta = Inf, control = list(l.chain = 3)), "'Sigma.beta'")
  expect_error(bayesGARCH(y, lambda = Inf, control = list(l.chain = 3)), "'lambda'")
  expect_error(bayesGARCH(y, delta = Inf, control = list(l.chain = 3)), "'delta' must be finite")
  expect_error(bayesGARCH(y, mu.alpha = c(0, Inf), control = list(l.chain = 3)), "'mu.alpha' must be finite")
  expect_error(bayesGARCH(y, Sigma.alpha = matrix(c(1000, 0, 0, Inf), 2),
                          control = list(l.chain = 3)), "'Sigma.alpha' must be finite")
})

test_that("formSmpl rejects malformed chain lists and non-integer settings", {
  M <- matrix(seq_len(40), 10, 4, dimnames = list(NULL, c("alpha0", "alpha1", "beta", "nu")))
  expect_error(formSmpl(list(), l.bi = 1), "empty")
  expect_error(formSmpl(list(M, M[, 1:3]), l.bi = 2), "same number of columns")
  expect_error(formSmpl(list(M, M[1:5, ]), l.bi = 2), "same number of rows")
  expect_error(formSmpl(list(M, unname(M)), l.bi = 2), "same column names")
  expect_error(formSmpl(list(M, "a"), l.bi = 2), "must be a matrix")
  expect_error(formSmpl(M, l.bi = 2.5), "whole number")
  expect_error(formSmpl(M, l.bi = Inf), "whole number")
  expect_error(formSmpl(M, l.bi = 2, batch.size = Inf), "whole number")
  expect_error(formSmpl(M, l.bi = 2, batch.size = 1.5), "whole number")
})

test_that("addPriorConditions is validated identically at the start and at proposals", {
  ok  <- function(psi) psi[2] + psi[3] < 1
  ## a function returning a truthy number is legitimate and was rejected up front
  ## while being accepted at proposal time
  num <- function(psi) as.numeric(psi[2] + psi[3] < 1)
  expect_silent(bayesGARCH(y, control = list(l.chain = 5, refresh = 0,
                                             addPriorConditions = num)))
  m <- bayesGARCH(y, control = list(l.chain = 20, refresh = 0,
                                    addPriorConditions = num))
  expect_true(all(as.matrix(m[[1]])[, "alpha1"] + as.matrix(m[[1]])[, "beta"] < 1))

  ## malformed results are refused with the same message wherever they occur
  expect_error(bayesGARCH(y, control = list(l.chain = 5, addPriorConditions = function(psi) NA)),
               "single non-missing")
  expect_error(bayesGARCH(y, control = list(l.chain = 5,
                 addPriorConditions = function(psi) c(TRUE, TRUE))), "single non-missing")
  expect_error(bayesGARCH(y, control = list(l.chain = 5,
                 addPriorConditions = function(psi) "yes")), "single non-missing")
  expect_error(bayesGARCH(y, control = list(l.chain = 5,
                 addPriorConditions = function(psi) stop("boom"))), "boom")

  ## a condition that only misbehaves at proposed states is reported the same way
  bad <- function(psi) if (psi[1] > 0 && psi[2] > 0.1000001) NA else TRUE
  expect_error(bayesGARCH(y, control = list(l.chain = 60, refresh = 0,
                 start.val = c(0.01, 0.1, 0.7, 10), addPriorConditions = bad)),
               "single non-missing")
})

test_that("priors too concentrated to invert are refused with a clear message", {
  expect_error(bayesGARCH(y, Sigma.beta = 1e-320, control = list(l.chain = 3)),
               "reciprocal is not representable")
  expect_error(bayesGARCH(y, lambda = 1e-320, control = list(l.chain = 3)),
               "reciprocal is not representable")
  expect_error(bayesGARCH(y, Sigma.alpha = diag(1e-320, 2), control = list(l.chain = 3)),
               "cannot be inverted|not representable")
})

test_that("covariance symmetry is judged to tolerance, not to the last bit", {
  ## the two off-diagonal entries are arithmetically equal but not bit-identical
  S <- matrix(c(1000, 0.1 + 0.2, 0.3, 1000), 2)
  expect_false(S[1, 2] == S[2, 1])
  expect_silent(bayesGARCH(y, Sigma.alpha = S, control = list(l.chain = 3, refresh = 0)))
  ## genuinely asymmetric matrices are still refused
  expect_error(bayesGARCH(y, Sigma.alpha = matrix(c(1000, 2, 3, 1000), 2),
                          control = list(l.chain = 3)), "not symmetric")
})
