\name{formSmpl}
\alias{formSmpl}
\title{Form the Posterior Sample}
\description{%
  Form the joint posterior sampler from the MCMC output.
}
\usage{%
  formSmpl(MCMC, l.bi = 0, batch.size = 1)
}
\arguments{%
  \item{MCMC}{object of the class \code{mcmc.list} (\R package
    \pkg{coda}) or a list of matrices or a matrix.}
  \item{l.bi}{length of the \emph{burn-in} phase.}
  \item{batch.size}{batching size used to diminish the autocorrelation within the chains.}
}
\value{%
  The joint posterior sample as an \code{mcmc} object (\R package \pkg{coda}).
}
\note{%
  Please cite the package in publications. Use \code{citation("bayesGARCH")}.
}
\seealso{%
  \code{\link{bayesGARCH}} for the Bayesian estimation of the GARCH(1,1)
  model with Student-t innovations.
}
\examples{
  ## !!! INCREASE THE NUMBER OF MCMC ITERATIONS !!!

  ## LOAD DATA SET
  data(dem2gbp)
  y <- dem2gbp[1:750]

  ## RUN THE ESTIMATION
  MCMC <- bayesGARCH(y, control = list(n.chain = 2, l.chain = 100))

  ## FORM THE SAMPLE FROM THE MCMC OUTPUT
  smpl <- formSmpl(MCMC, l.bi = 50, batch.size = 2)

  ## POSTERIOR STATISTICS
  summary(smpl)
}
\keyword{misc}
