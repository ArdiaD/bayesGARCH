## fn.dmvnorm2() and fn.rmvnorm2() replaced mvtnorm::dmvnorm() and
## mvtnorm::rmvnorm(method = "eigen") in the alpha step. They are meant to be
## exact substitutes, not approximations: same value, and for the draw, the same
## value from the same seed. If that ever stops holding, the sampler's draws
## change silently, so the equivalence is pinned here.

test_that("fn.dmvnorm2 reproduces mvtnorm::dmvnorm exactly", {
  skip_if_not_installed("mvtnorm")
  fn.dmvnorm2 <- bayesGARCH:::fn.dmvnorm2

  set.seed(20260810)
  for (k in 1:200) {
    A <- matrix(rnorm(4), 2)
    sigma <- solve(crossprod(A) + diag(2) * runif(1, 0.1, 2))
    mu <- rnorm(2)
    x  <- rnorm(2)

    target <- unname(mvtnorm::dmvnorm(x, mean = mu, sigma = sigma, log = TRUE))
    expect_identical(fn.dmvnorm2(x, mean = mu, sigma = sigma), target)
    ## the alpha step passes a 2x1 column and a 2x1 mean
    expect_identical(fn.dmvnorm2(matrix(x, 2, 1), mean = matrix(mu, 2, 1), sigma = sigma),
                     structure(target, dim = NULL))
  }
})

test_that("fn.rmvnorm2 reproduces mvtnorm::rmvnorm(method = 'eigen') draw for draw", {
  skip_if_not_installed("mvtnorm")
  fn.rmvnorm2 <- bayesGARCH:::fn.rmvnorm2

  for (k in 1:200) {
    set.seed(k)
    A <- matrix(rnorm(4), 2)
    sigma <- solve(crossprod(A) + diag(2) * runif(1, 0.1, 2))
    mu <- rnorm(2)

    set.seed(1000 + k)
    target <- as.numeric(mvtnorm::rmvnorm(1, mean = mu, sigma = sigma))
    set.seed(1000 + k)
    expect_identical(as.numeric(fn.rmvnorm2(mean = mu, sigma = sigma)), target)
  }
})

## fn.Dd.D() builds the posterior precision from the whitened regressors
## Z := X/sqrt(h). Going through crossprod(Z) rather than t(X/h) %*% X makes the
## result symmetric by construction rather than by luck, which matters because
## the matrix is handed straight to chol() and to eigen(symmetric = TRUE), both
## of which read a single triangle and would silently ignore any asymmetry.

test_that("fn.Dd.D returns an exactly symmetric D and the right cross-products", {
  fn.Dd.D <- bayesGARCH:::fn.Dd.D

  set.seed(20260810)
  for (k in 1:100) {
    n <- 500
    X  <- cbind(runif(n, 0.5, 7), runif(n, 0, 13))
    h  <- runif(n, 5e-5, 8)
    yv <- runif(n, 0, 40)
    iv <- solve(1000 * diag(1, 2))
    r  <- fn.Dd.D(yv, X, h, c(0, 0), iv)

    expect_identical(r$D[1, 2], r$D[2, 1])
    expect_equal(solve(r$D), t(X / h) %*% X + iv, tolerance = 1e-12)
    expect_equal(as.numeric(r$Dd), as.numeric(r$D %*% (t(X / h) %*% yv)),
                 tolerance = 1e-10)
  }

  ## the beta step passes plain vectors and a scalar prior precision
  n <- 500
  W <- rnorm(n); h <- runif(n, 0.1, 5); z <- rnorm(n)
  r <- fn.Dd.D(z, W, h, 0, 1/1000)
  expect_equal(dim(r$D), c(1L, 1L))
  expect_equal(drop(solve(r$D)), sum(W * W / h) + 1/1000, tolerance = 1e-12)
})
