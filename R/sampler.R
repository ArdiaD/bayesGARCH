## sampler.R
## David Ardia

## bayesGARCH
## __input__
## y : time series
## mu.alpha = c(0,0) : prior mean for alpha
## Sigma.alpha = 1000 * diag(1,2) : prior covariance matrix for alpha
## mu.beta = 0 : prior mean for beta
## Sigma.beta = 1000 : prior variance for beta
## lambda = 0.01 : rate of the translated Exponential prior on nu
## delta = 2 : lower bound of that prior; the prior mean of nu is delta + 1/lambda
## control = list() : a list with
##           n.chain : number of chains; defaults to nrow(start.val)
##           l.chain = 10000 : length of each chain
##           start.val = c(0.01, 0.1, 0.7, 100) : starting values, or one row per chain
##           addPriorConditions : a function of psi returning TRUE if psi is admissible
##           digits = 4 : digits for the output
##           refresh = 10 : refreshing frequency
## __output__
## r : an object of class mcmc.list of length n.chain containing the MCMC outputs
"bayesGARCH" <- function(y, mu.alpha = c(0,0), Sigma.alpha = 1000 * diag(1,2),
                         mu.beta = 0, Sigma.beta = 1000,
                         lambda = 0.01, delta = 2, control = list()){
  if (missing(y))
    stop ("'y' is missing")
  ## a 'ts' or 'zoo' series is not a vector in the is.vector() sense, but it is
  ## the natural input for this function
  if (!is.numeric(y) || !is.null(dim(y)))
    stop ("'y' must be a numeric vector")
  y <- as.numeric(y)
  if (length(y) < 2)
    stop ("'y' must be a longer time series")
  if (any(is.na(y)))
    stop ("'y' contains 'NA' values")
  if (any(!is.finite(y)))
    stop ("'y' contains non-finite values")
  if (!is.vector(mu.alpha) || !is.numeric(mu.alpha))
    stop ("'mu.alpha' must be a numeric vector")
  if (length(mu.alpha) != 2) 
    stop ("'mu.alpha' is not of appropriate size")
  if (any(!is.finite(mu.alpha)))
    stop ("'mu.alpha' must be finite")
  if (!is.matrix(Sigma.alpha) || !is.numeric(Sigma.alpha))
    stop ("'Sigma.alpha' must be a numeric matrix")
  if (!all(dim(Sigma.alpha) == c(2, 2))) 
    stop ("'Sigma.alpha' is not of appropriate size")
  if (any(!is.finite(Sigma.alpha)))
    stop ("'Sigma.alpha' must be finite")
  if (Sigma.alpha[1,2] != Sigma.alpha[2,1])
    stop ("'Sigma.alpha' is not symmetric")
  if (any(eigen(Sigma.alpha, symmetric = TRUE, only.values = TRUE)$values <= 0))
    stop ("'Sigma.alpha' is not positive definite")
  if (!is.vector(mu.beta) || !is.numeric(mu.beta) || length(mu.beta) != 1)
    stop ("'mu.beta' must be a scalar")
  if (!is.finite(mu.beta))
    stop ("'mu.beta' must be finite")
  if (!is.vector(Sigma.beta) || !is.numeric(Sigma.beta) || length(Sigma.beta) != 1)
    stop ("'Sigma.beta' must be a scalar")
  if (!is.finite(Sigma.beta) || Sigma.beta <= 0)
    stop ("'Sigma.beta' must be a finite positive value")
  if (!is.vector(lambda) || !is.numeric(lambda) || length(lambda) != 1)
    stop ("'lambda' must be a scalar")
  if (!is.finite(lambda) || lambda <= 0)
    stop ("'lambda' must be a finite positive value")
  if (!is.vector(delta) || !is.numeric(delta) || length(delta) != 1)
    stop ("'delta' must be a scalar")
  if (!is.finite(delta))
    stop ("'delta' must be finite")
  if (delta < 2)
    stop ("'delta' cannot be lower than '2'")
  if (!is.list(control))
    stop ("'control' must be a list")
  if (!is.null(control$addPriorConditions))
    if (!is.function(control$addPriorConditions))
      stop ("'control$addPriorConditions' must be NULL or a function")
  
  fn.bayesGARCH(y, mu.alpha, solve(Sigma.alpha), mu.beta, 1/Sigma.beta, lambda, delta, control)
}

## fn.checkCount
## a control parameter that counts something: a single finite whole number
"fn.checkCount" <- function(x, name, min){
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x != round(x) || x < min)
    stop ("'", name, "' must be a single finite whole number greater or equal than '", min, "'")

  invisible(TRUE)
}

## fn.checkStartVal
## Every chain must start inside the support of the posterior it is meant to
## explore. The Metropolis-Hastings steps only ever accept a move into the
## admissible region, so a chain started outside it can never enter: the
## starting value is returned unchanged at every iteration, and the whole sample
## is made of draws that the target gives zero mass. Rather than let that pass
## silently, the starting values are checked up front.
"fn.checkStartVal" <- function(start.val, d.nu, addPriorConditions){
  nms <- c("alpha0", "alpha1", "beta", "nu")
  for (i in 1:nrow(start.val)) {
    psi <- start.val[i,]
    lab <- paste0("'control$start.val' row ", i)
    if (any(!is.finite(psi)))
      stop (lab, " contains non-finite values")
    if (any(psi[1:3] <= 0))
      stop (lab, ": ", paste(nms[1:3][psi[1:3] <= 0], collapse = ", "),
            " must be positive")
    if (psi[4] <= d.nu)
      stop (lab, ": 'nu' must be greater than 'delta' (", format(d.nu), ")")
    ok <- try(addPriorConditions(psi), silent = TRUE)
    if (inherits(ok, "try-error"))
      stop (lab, ": 'control$addPriorConditions' failed when applied to it")
    if (!isTRUE(ok))
      stop (lab, " does not satisfy 'control$addPriorConditions'")
  }

  invisible(TRUE)
}

## sub function fn.bayesGARCH
"fn.bayesGARCH" <- function(y, mu.alpha, iv.alpha, mu.beta, iv.beta, c.nu, d.nu, control){

  ## the default starting value for nu is the prior mean delta + 1/lambda, so that
  ## it lies in the prior support whatever 'delta' is. A fixed value would sit
  ## outside it for the Normal special case, where delta = 500.
  con <- list(n.chain = NULL, l.chain = 10000, addPriorConditions = NULL, digits = 4, refresh = 10,
              start.val = matrix(c(0.01, 0.1, 0.7, d.nu + 1/c.nu), 1, 4,
                dimnames = list("chain1", c("alpha0", "alpha1", "beta", "nu"))),
              hypers = list(mu.alpha = mu.alpha, iv.alpha = iv.alpha,
                mu.beta = mu.beta, iv.beta = iv.beta, c.nu = c.nu, d.nu = d.nu))

  ## 'hypers' is internal: the hyper-parameters are validated as arguments of
  ## bayesGARCH() and must not be replaceable through 'control'.
  known <- c("n.chain", "l.chain", "start.val", "addPriorConditions", "refresh", "digits")
  if (is.null(names(control)) && length(control) > 0L)
    stop ("'control' must be a named list")
  unknown <- setdiff(names(control), known)
  if (length(unknown) > 0L)
    stop ("unknown component(s) in 'control': ", paste(sQuote(unknown), collapse = ", "),
          ". Valid components are ", paste(sQuote(known), collapse = ", "), ".")
  con[names(control)] <- control

  fn.checkCount(con$l.chain, "control$l.chain", 2)
  fn.checkCount(con$refresh, "control$refresh", 0)
  fn.checkCount(con$digits,  "control$digits",  0)
  if (is.vector(con$start.val)) 
    con$start.val <- matrix(con$start.val, nrow = 1, byrow = TRUE)
  if (!is.numeric(con$start.val))
    stop ("'control$start.val' must be numeric")
  con$start.val <- as.matrix(con$start.val)
  if (ncol(con$start.val) != 4)
    stop ("'control$start.val' must have four columns ('alpha0', 'alpha1', 'beta', 'nu')")
  if (nrow(con$start.val) < 1)
    stop ("'control$start.val' must have at least one row")
  if (is.null(con$n.chain)) 
    con$n.chain <- nrow(con$start.val)
  ## a single row of starting values is replicated across the chains; several rows
  ## are used as they stand (one row per chain), which is what the documentation
  ## promises. Before 2.2.0 the second case silently recycled the whole matrix and
  ## then failed with an obscure 'dimnames' error.
  else {
    fn.checkCount(con$n.chain, "control$n.chain", 1)
    if (nrow(con$start.val) == 1)
      con$start.val <- matrix(con$start.val, con$n.chain, 4, byrow = TRUE)
    else if (nrow(con$start.val) != con$n.chain)
      stop ("'control$n.chain' does not match the number of rows of 'control$start.val'")
  }
  if (is.null(colnames(con$start.val))) 
    colnames(con$start.val) <- c("alpha0", "alpha1", "beta", "nu") 
  if (is.null(rownames(con$start.val))) 
    rownames(con$start.val) <- paste(rep("chain", con$n.chain), 1:con$n.chain, sep = "")
  if (is.null(con$addPriorConditions)) {con$addPriorConditions <- function(psi){TRUE}}
  fn.checkStartVal(con$start.val, d.nu, con$addPriorConditions)
  r <- list()
  for (i in 1:con$n.chain) {
    k <- NULL
    chain <- fn.initializeChain(con$l.chain, con$start.val[i,])
    for (j in 2:con$l.chain) {
      k <- j - 1
      chain[j, ] <- fn.block(y, chain[k, 1:2], chain[k,3], chain[k, 4], con$hypers$mu.alpha, con$hypers$iv.alpha, 
                             con$hypers$mu.beta, con$hypers$iv.beta, con$hypers$c.nu, 
                             con$hypers$d.nu, con$addPriorConditions)
      if (con$refresh > 0 & j%%con$refresh == 0) 
        message("chain: ", i, " iteration: ", j, " parameters: ", 
                paste(round(chain[j, ], con$digits), collapse = " "))
    }
    ## proportion of iterations at which each parameter moved. A high rate is not
    ## by itself evidence of good mixing: the accepted moves are typically small.
    attr(chain, "move.rates") <- colMeans(abs(diff(as.matrix(chain))) > 0)
    r[[i]] <- chain
  }
  names(r) <- paste(rep("chain", con$n.chain), 1:con$n.chain, sep = "")
  r <- do.call(coda::mcmc.list, r)
  names(r) <- paste(rep("chain", con$n.chain), 1:con$n.chain, sep = "")
  r
}

"fn.cvGARCH" <- function(y, theta){
  n <- length(y)
  .C('fnGarchC',
     n = as.integer(n),
     order = as.integer(c(1,2)),
     alpha = as.double(theta[1:2]),
     delta = as.double(0), 
     beta = as.double(theta[3]),
     h = vector('double', n+1),
     u = as.double(y),
     sim = vector('integer',n+1),
     PACKAGE = 'bayesGARCH')$h
}

## fn.initializeChain
## __input__
## l.chain : length of the chain
## start.val : starting values
## __output__
## mcmc matrix 
"fn.initializeChain" <- function (l.chain, start.val){
  r <- matrix(NA_real_, l.chain, length(start.val), dimnames = list(1:l.chain, names(start.val)))
  r[1,] <- start.val
  
  coda::mcmc(r, start = 1)
}

## fn.block
## __input__
## y : data
## alpha, beta, nu : parameters
## alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu : hyperparameters
## __output__
## block of parameters
"fn.block" <- function(y, alpha, beta, nu, alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu, addPriorConditions){
  
  ## the w and alpha steps both condition on h evaluated at the current
  ## (alpha, beta); compute it once
  h <- fn.cvGARCH(y, c(alpha, beta))[seq_along(y)]

  w.new <- fn.w.full(y, alpha, beta, nu,
                     alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu, h = h)	
  
  ## The prior constraint must enter EACH Metropolis step's acceptance test.
  ## Applying it once to the whole block does not leave p(psi|y) I[psi in A]
  ## invariant: the block is a composition of pi-reversible kernels, which is
  ## pi-invariant but not itself pi-reversible.
  alpha.new <- fn.alpha.full(y, alpha, beta, nu, w.new,
                             alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu, addPriorConditions,
                             h = h)

  beta.new <- fn.beta.full(y, alpha.new, beta, nu, w.new,
                           alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu, addPriorConditions)

  nu.new <- fn.nu.full(y, alpha.new, beta.new, nu, w.new,
                       alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu, addPriorConditions)

  c(alpha.new, beta.new, nu.new)
}

## fn.alpha.full
## __input__
## y : data
## alpha.old, beta, nu, w : parameters
## alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu : hyperparameters
## __output__
## MH draw for alpha
"fn.alpha.full" <- function(y, alpha.old, beta, nu, w,
                            alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu,
                            addPriorConditions = function(psi) TRUE, h = NULL){
  r <- alpha.old
  n <- length(y)
  if (is.null(h))
    h <- fn.cvGARCH(y, c(alpha.old, beta))[1:n]
  
  post.old <- fn.post.garch(y, w*h, alpha.old, beta, nu,
                            alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu)

  tau <- 2 * h^2
  filter.alpha <- fn.filterAlpha(y, beta)

  Dd.D <- fn.Dd.D(y^2/w, filter.alpha, tau,
                  alpha0, iv.alpha0)
  
  alpha.new <- fn.rmvnorm2(mean = Dd.D$Dd, sigma = Dd.D$D)

  if (all(alpha.new>0)){
    h.new <- fn.cvGARCH(y, c(alpha.new, beta))[1:n]

    tau.new <- 2 * h.new^2
		
    post.new <- fn.post.garch(y, w*h.new, alpha.new, beta, nu,
                              alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu)

    prop.new <- fn.dmvnorm2(alpha.new, mean = Dd.D$Dd, sigma = Dd.D$D)
		
    Dd.D <- fn.Dd.D(y^2/w, filter.alpha, tau.new,
                    alpha0, iv.alpha0)
		
    prop.old <- fn.dmvnorm2(alpha.old, mean = Dd.D$Dd, sigma = Dd.D$D)

    ratio <- post.new - post.old + prop.old - prop.new

    r.u <- runif(1)
    if (isTRUE(r.u < min(1, exp(ratio))) && addPriorConditions(c(alpha.new, beta, nu)))
      r <- alpha.new
  }
  
  r
}

## fn.beta.full
## __input__
## y : data set
## alpha, beta.old, nu, w : parameters
## alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu : hyperparameters
## __output__
## MH draw for beta
"fn.beta.full" <- function(y, alpha, beta.old, nu, w,
                           alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu,
                           addPriorConditions = function(psi) TRUE){
  n <- length(y)
  r <- beta.old
  h <- fn.cvGARCH(y, c(alpha, beta.old))[1:n]

  post.old <- fn.post.garch(y, w*h, alpha, beta.old, nu,
                            alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu)
  tau <- 2 * h^2
  beta.nls <- beta.old

  w1 <- fn.filterW(y, alpha, beta.nls)
  W <- fn.W(y^2-w1, -beta.nls)
  z <- w1 + W * beta.nls
  z <- z + y^2 * (1/w - 1)
  
  Dd.D <- fn.Dd.D(z, W, tau, beta0, iv.beta0)

  beta.new <- rnorm(n = 1, mean = Dd.D$Dd, sd = sqrt(Dd.D$D))
	
  if (beta.new>0){
    h.new <- fn.cvGARCH(y, c(alpha, beta.new))[1:n]
    tau.new <- 2 * h.new^2
    
    post.new <- fn.post.garch(y, w*h.new, alpha, beta.new, nu, 
                              alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu)
				
    prop.new <- dnorm(beta.new, mean = Dd.D$Dd, sd = sqrt(Dd.D$D), log = TRUE)
		
    w1 <- fn.filterW(y, alpha, beta.new)
    W <- fn.W(y^2-w1, -beta.new)
    z <- w1 + W * beta.new
    z <- z + y^2 * (1/w - 1)
    
    Dd.D <- fn.Dd.D(z, W, tau.new, beta0, iv.beta0)
		
    prop.old <- dnorm(beta.old, mean = Dd.D$Dd, sd = sqrt(Dd.D$D), log = TRUE)

    ratio <- post.new - post.old + prop.old - prop.new

    r.u <- runif(1)
    if (isTRUE(r.u < min(1, exp(ratio))) && addPriorConditions(c(alpha, beta.new, nu)))
      r <- beta.new  
  }
  
  r  	
}

## fn.w.full
## __input__
## u : data set
## alpha, beta, nu : parameters
## alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu : hyperparameters
## __output__
## draw for w
"fn.w.full" <- function(u, alpha, beta, nu,
                        alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu, h = NULL){
  n <- length(u)
  if (is.null(h))
    h <- fn.cvGARCH(u, c(alpha, beta))[1:n]

  ## the shape is the same for every t; rgamma() recycles it
  r <- rgamma(n = n, shape = 0.5 * (nu+1),
              rate = 0.5 * ((u^2/h) + (nu-2)) )
  
  1/r
}

## fn.nu.full
## __input__
## u : data set
## alpha, beta, nu.old, w : parameters
## alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu : hyperparameters
## alpha.min, alpha.max, k.max : simulation parameters
## __output__
## simulated nu
"fn.nu.full" <- function(u, alpha, beta, nu.old, w,
                         alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu,
                         addPriorConditions = function(psi) TRUE,
                         alpha.min = 0.00001, alpha.max = 500, k.max = 50000){
  n <- length(u)
  r <- nu.old
  alpha.garch <- alpha          ## `alpha` below is reused for the exponential rate
		
  phi <- 0.5 * sum(log(w) + 1/w) + c.nu
  
  optim.alpha <- uniroot(fn.neg.alpha,
                         interval = c(alpha.min, alpha.max),
                         tol = .Machine$double.eps,
                         maxiter = 1000,
                         n = n,
                         d = d.nu,
                         phi = phi)
  
  alpha <- optim.alpha$root
  if ( (alpha>alpha.min) && (alpha<alpha.max) ){
    k <- 0
    while (k<k.max){
      k <- k+1
      nu.new <- rexp(n = 1, rate = alpha)
      nu.new <- nu.new + d.nu
       
      p <- fn.accept.proba(n, alpha, d.nu, nu.new, phi)
    
      r.u <- runif(1)
      if(isTRUE(r.u < p)){
        ## nu.new is now an exact draw from the full conditional, so imposing the
        ## prior constraint is an independence Metropolis-Hastings step whose
        ## ratio reduces to the indicator: keep the draw if it is admissible and
        ## the current value otherwise. Testing the constraint inside the
        ## accept-reject loop instead would spin for all 'k.max' proposals
        ## whenever no value of nu can satisfy it, which is the usual case since
        ## conditions such as alpha1 + beta < 1 do not involve nu at all.
        if (addPriorConditions(c(alpha.garch, beta, nu.new))) r <- nu.new
        k = k.max
      }
    }
  }
  
  r
}

## fn.post.garch
## __input__
## u : data set
## h : conditional variance process
## alpha, beta, nu : parameters
## alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu : parameters
## __output__
## posterior density
"fn.post.garch" <- function(u, h, alpha, beta, nu,
                            alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu){

  r <- 0
  r <- -0.5 * sum(log(h))
  r <- r - 0.5 * sum(u^2/h)
  r <- r - 0.5 * t(alpha-alpha0) %*%
    iv.alpha0 %*% (alpha-alpha0)
  r <- r - 0.5 * (beta-beta0)^2 * iv.beta0
  r <- r - c.nu * nu
  
  r
}

## fn.neg.alpha
"fn.neg.alpha" <- function(alpha, n, d, phi){
  r <- 0
  r <- log( (1+alpha*(d-2))/(2*alpha) )
  r <- r + (1+alpha*d)/(1+alpha*(d-2))
  r <- r - digamma( (1+alpha*d)/(2*alpha) )
  r <- 0.5 * n * r + alpha - phi
  
  r
}

## fn.accept.proba
"fn.accept.proba" <- function(n, alpha, d, nu, phi){
  u <- (1+alpha*d) / (2*alpha)
  v <- nu/2
  w <- (1+alpha*(d-2)) / (2*alpha)
  
  r.ln <- 0
  r.ln <- n * (lgamma(u) - lgamma(v))
  r.ln <- r.ln + n * v * log((nu-2)/2)
  r.ln <- r.ln - n * u * log(w)
  r.ln <- r.ln + (nu - d) * (alpha - phi) + (phi/alpha) - 1
  
  exp(r.ln)
}

## fn.filterAlpha
"fn.filterAlpha" <- function(u, beta){
  n <- length(u)	
  filter.alpha <- .C('fnFilterAlphaC', 
                     n = as.integer(n),
                     u = as.double(u),
                     beta = as.double(beta),
                     vstar = vector('double', n),
                     lstar = vector('double', n),
                     PACKAGE = 'bayesGARCH')
  
  as.matrix(cbind(filter.alpha$lstar, 
                  c(0, filter.alpha$vstar[1:(n-1)])))
}

## fn.filterW
"fn.filterW" <- function(u, alpha, beta){
  n <- length(u)
  .C('fnFilterWC', 
     n = as.integer(n),
     u = as.double(u),
     alpha = as.double(alpha),
     beta = as.double(beta),
     w = vector('double', n),
     PACKAGE = 'bayesGARCH')$w
}

## fn.W
"fn.W" <- function(u, theta){
  n <- length(u)
  .C('fnQDiffC',
     n = as.integer(n),
     u = as.double(u),
     theta = as.double(theta),
     w = vector('double', n),
     w0 = as.double(0),
     PACKAGE = 'bayesGARCH')$w
}

## fn.dmvnorm2 / fn.rmvnorm2
## Log-density of, and a draw from, a bivariate Normal. These reproduce
## mvtnorm::dmvnorm() and mvtnorm::rmvnorm(method = "eigen") operation for
## operation, and therefore draw for draw from a given seed, without the
## argument checking and dispatch that dominate their cost at p = 2.
"fn.dmvnorm2" <- function(x, mean, sigma){
  dec <- chol(sigma)
  tmp <- backsolve(dec, x - mean, transpose = TRUE)

  -sum(log(diag(dec))) - log(2*pi) - 0.5 * sum(tmp^2)
}

"fn.rmvnorm2" <- function(mean, sigma){
  ev <- eigen(sigma, symmetric = TRUE)
  R <- t( ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))) )

  drop( matrix(rnorm(2), nrow = 1) %*% R ) + mean
}

## fn.Dd.D
## compute D and Dd
## The weighted cross-products are formed from the whitened regressors
## Z := X / sqrt(h) rather than from X/h. Both give X' diag(1/h) X, but
## crossprod(Z) goes through a symmetric rank-k update, so the result is
## symmetric by construction and about five times closer to the exactly
## rounded value; crossprod(Z, y/s) is likewise more accurate than
## t(X/h) %*% y when the regressors change sign, as they do in the beta
## step. It is also the cheaper of the two.
##
## The posterior precision is symmetric positive definite -- the prior
## contributes iv.b0 -- so it is inverted through its Cholesky factor. Beyond
## being faster than a general LU solve, chol2inv() returns an exactly
## symmetric D, which solve() does not: it did so for fewer than half of the
## matrices arising here. That matters because D is passed on to chol(), which
## reads its upper triangle, and to eigen(symmetric = TRUE), which reads its
## lower one. With an asymmetric D the covariance used to draw the candidate
## and the covariance used to evaluate its density are not quite the same
## matrix, and the Metropolis-Hastings ratio is correspondingly inconsistent.
## The discrepancy is of the order of the rounding error and of no practical
## consequence, but it costs nothing to remove.
"fn.Dd.D" <- function(y, X, h, b0, iv.b0){
  s <- sqrt(h)
  Z <- X/s
  D <- chol2inv(chol(crossprod(Z) + iv.b0))
  d <- crossprod(Z, y/s) + iv.b0 %*% b0

  list(D = D, Dd = D %*% d)
}


