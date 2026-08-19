## LIBRARY garchLib
## functions.bayes
## David Ardia

## formSmpl
## __input__
## MCMC : and object of the class MCMC.list or a list of matrices containing the MCMC output
## l.bi = 0 : the lenth of the burn in phase
## batch.size : size of the batching smpl
## __output__
## r : a matrix of the posterior sample
## formSmpl
## __input__
## MCMC : and object of the class MCMC.list or a list of matrices containing the MCMC output
## l.bi = 0 : the lenth of the burn in phase
## batch.size = 1 : size of the batching smpl
## __output__
## r : a matrix of the posterior sample (object of class mcmc)
"formSmpl" <- function(MCMC, l.bi = 0, batch.size = 1){
  if (missing(MCMC))
    stop ("'MCMC' is missing")
  if (!is.list(MCMC)){
    if (is.matrix(MCMC))
      MCMC <- list(MCMC)
    else
      stop ("'MCMC' is neither a list or a matrix")
  }
  if (!is.vector(l.bi) || !is.numeric(l.bi) || length(l.bi) != 1)
    stop ("'l.bi' must be a scalar")
  if (!is.finite(l.bi) || l.bi != round(l.bi) || l.bi < 0)
    stop ("'l.bi' must be a finite non-negative whole number")
  if (!is.vector(batch.size) || !is.numeric(batch.size) || length(batch.size) != 1)
    stop ("'batch.size' must be a scalar")
  if (!is.finite(batch.size) || batch.size != round(batch.size) || batch.size < 1)
    stop ("'batch.size' must be a finite whole number larger or equal than 1")
  r <- NULL
  n.chain <- length(MCMC)
  if (n.chain < 1)
    stop ("'MCMC' is empty")
  ## every chain must be a matrix-like object of the same width, otherwise the
  ## chains are bound together by position and the result is silently wrong
  if (!all(vapply(MCMC, function(z) is.matrix(z) || is.data.frame(z), NA)))
    stop ("every element of 'MCMC' must be a matrix")
  ncols <- vapply(MCMC, ncol, 1L)
  if (any(ncols != ncols[1]))
    stop ("all chains in 'MCMC' must have the same number of columns")
  cn <- lapply(MCMC, colnames)
  if (!all(vapply(cn, function(z) identical(z, cn[[1]]), NA)))
    stop ("all chains in 'MCMC' must have the same column names")
  l.chain <- nrow(MCMC[[1]])
  if (any(vapply(MCMC, nrow, 1L) != l.chain))
    stop ("all chains in 'MCMC' must have the same number of rows")
  if (l.bi >= l.chain)
    stop ("'l.bi >= l.chain'")
  theM <- ((l.bi + 1):l.chain)
  seqk <- seq(from = 1, to = length(theM), by = batch.size)
  theM <- theM[seqk]
  for (i in 1:n.chain)
    r <- rbind(r, MCMC[[i]][theM, ])
  message("n.chain: ", n.chain, "  l.chain: ", l.chain, "  l.bi: ", l.bi,
          "  batch.size: ", batch.size, "  smpl size: ", nrow(r))
  dimnames(r) <- list(1:nrow(r), colnames(MCMC[[1]]))

  coda::mcmc(r, start = 1)
}
