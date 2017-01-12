\name{bayesGARCH}
\alias{bayesGARCH}
\title{Bayesian Estimation of the GARCH(1,1) Model with Student-t Innovations}
\concept{Bayesian GARCH}
\description{%
  Performs the Bayesian estimation of the GARCH(1,1) model with
  Student-t innovations.
}
\usage{%
  bayesGARCH(y, mu.alpha = c(0,0), Sigma.alpha = 1000 * diag(1,2), 
             mu.beta = 0, Sigma.beta = 1000,
             lambda = 0.01, delta = 2, control = list())
}
\arguments{%
  \item{y}{vector of observations of size \eqn{T}. \code{NA} values are not allowed.}
  \item{mu.alpha}{hyper-parameter \eqn{\mu_\alpha}{mu_alpha} (prior mean)
    for the truncated Normal prior on parameter
    \eqn{\alpha := (\alpha_0 \ \alpha_1)'}{alpha:=(alpha0
      alpha1)'}. Default: a \eqn{2 \times 1}{2x1} vector of zeros.}
  \item{Sigma.alpha}{hyper-parameter \eqn{\Sigma_\alpha}{Sigma_alpha}
    (prior covariance matrix) for
    the truncated Normal prior on parameter
    \eqn{\alpha}{alpha}. Default: a \eqn{2 \times 2}{2x2} diagonal matrix whose
    variances are set to 1'000, i.e., a diffuse prior. Note that the matrix
    must be symmetric positive definite.}
  \item{mu.beta}{hyper-parameter \eqn{\mu_\beta}{mu_beta} (prior mean) for the truncated Normal prior on
    parameter \eqn{\beta}{beta}. Default: zero.}
  \item{Sigma.beta}{hyper-parameter \eqn{\Sigma_\beta > 0}{Sigma_beta>0} (prior variance) for the truncated
    Normal prior on parameter \eqn{\beta}{beta}. Default: 1'000, i.e.,
    a diffuse prior.} 
  \item{lambda}{hyper-parameter \eqn{\lambda > 0}{lambda>0} for the translated
    Exponential distribution on parameter \eqn{\nu}{nu}. Default: 0.01.}
  \item{delta}{hyper-parameter \eqn{\delta \ge 2}{delta>=2} for the translated
    Exponential distribution on parameter \eqn{\nu}{nu}. Default: 2 (to
    ensure the existence of the conditional variance).}
  \item{control}{list of control parameters (See *Details*).}
}
\details{%
  The function \code{bayesGARCH} performs the Bayesian estimation of the
  GARCH(1,1) model with Student-t innovations. The underlying algorithm is based on Nakatsuma
  (1998, 2000) for generating the parameters of the GARCH(1,1) scedastic
  function \eqn{\alpha := (\alpha_0 \ \alpha_1)'}{alpha:=(alpha0 alpha1)'} and \eqn{\beta}{beta} and on
  Geweke (1993) and Deschamps (2006) for the generating the degrees of freedom
  parameter \eqn{\nu}{nu}. Further details and examples can be found in Ardia (2008) and 
  Ardia and Hoogerheide (2010). See also the package vignette by typing \code{vignette("bayesGARCH")}. Finally, we refer to 
  Ardia (2009) for an extension of the algorithm to Markov-switching GARCH models.

  The \code{control} argument is a list that can supply any of
  the following components:

  \describe{%
    \item{\code{n.chain}}{number of MCMC chain(s) to be
      generated. Default: \code{n.chain=1}.}
    \item{\code{l.chain}}{length of each MCMC chain. Default: \code{l.chain=10000}.}    
    \item{\code{start.val}}{vector of starting values of
      chain(s). Default: \code{start.val=c(0.01,0.1,0.7,20)}. A matrix of
      size \eqn{n \times 4}{nx4}
      containing starting values in rows can also be provided. This will generate \eqn{n} chains starting at the
      different row values.}
    \item{\code{addPriorConditions}}{function which allows the user to add constraints on the model parameters. 
      Default: \code{NULL}, i.e. not additional constraints are imposed (see below).}
    \item{\code{refresh}}{frequency of reports. Default: \code{refresh=10} iterations.}
    \item{\code{digits}}{number of printed digits in the
      reports. Default: \code{digits=4}.}
  }
}
\value{%
  A list of class \code{mcmc.list} (\R package \pkg{coda}).
}
\note{%
  The GARCH(1,1) model with Student-t innovations may be written as follows:
  \deqn{
    y_t = \epsilon_t (\varrho h_t)^{1/2}
  }{
    y(t) = e(t)*(varrho * h(t))^(1/2)  
  }
  for \eqn{t=1,\ldots,T}{t=1,...,T}, where the conditional variance equation is defined as:
  \deqn{
    h_t := \alpha_0 + \alpha_1 y_{t-1}^2 + \beta h_{t-1}
  }{
    h(t) := alpha0 + alpha1 * y(t-1)^2 + beta * h(t-1)
  }
  where \eqn{\alpha_0 > 0, \alpha_1 \ge 0, \beta \ge 0}{alpha0>0,alpha1,beta>=0} to ensure a
  positive conditional variance. We set the initial variance to
  \eqn{h_0 := 0}{h(0):=0} for convenience. The parameter \eqn{\varrho :=
    (\nu-2)/\nu}{varrho:=(nu-2)/nu} is a scaling factor which ensures
  the conditional variance of \eqn{y_t}{y(t)} to be
  \eqn{h_t}{h(t)}. Finally, \eqn{\epsilon_t}{e(t)}
  follows a Student-t distribution with \eqn{\nu}{nu} degrees of freedom.

  The prior distributions on \eqn{\alpha}{alpha} is a bivariate truncated
  Normal distribution:
  \deqn{
    p(\alpha) \propto N_2(\alpha \mid \mu_\alpha, \Sigma_\alpha)
    I_{[\alpha > 0]}
  }{
    p(alpha) prop N2(alpha | mu_alpha, Sigma_alpha) I[alpha>0]
  }
  where \eqn{\mu_\alpha}{mu_alpha} is the prior mean vector, \eqn{\Sigma_\alpha}{Sigma_alpha} is
  the prior covariance matrix and \eqn{I_{[\bullet]}}{I[alpha>0]} is the indicator function.
  
  The prior distribution on \eqn{\beta}{beta} is a univariate truncated Normal
  distribution:
  \deqn{
    p(\beta) \propto N(\beta \mid \mu_\beta, \Sigma_\beta) I_{[\beta > 0]}
  }{
    p(theta) prop N(beta | mu_beta, Sigma_beta) I[beta>0]
  }
  where \eqn{\mu_\beta}{mu_beta} is the prior mean and \eqn{\Sigma_\beta}{Sigma_beta} is the
  prior variance.
  
  The prior distribution on \eqn{\nu}{nu} is a translated Exponential
  distribution:
  \deqn{
    p(\nu) = \lambda \exp [ - \lambda (\nu - \delta) ] I_{[\nu > \delta]}
  }{
    p(nu) = lambda * exp(-lambda(nu-delta)) I[nu>delta]
  }
  where \eqn{\lambda > 0}{lambda>0} and \eqn{\delta \ge 2}{delta>=2}. The prior mean for
  \eqn{\nu}{nu} is \eqn{\delta + 1/\lambda}{delta + 1/lambda}.
  
  The joint prior on parameter \eqn{\psi := (\alpha,
    \beta, \nu)}{psi:=(alpha,beta,nu)} is obtained by assuming prior independence:
  \deqn{
    p(\psi) = p(\alpha) p(\beta) p(\nu).
  }{
    p(psi) = p(alpha) * p(beta) * p(nu).
  }

  The default hyperparameters \eqn{\mu_\alpha}{mu_alpha}, \eqn{\Sigma_\alpha}{Sigma_alpha},
  \eqn{\mu_\beta}{mu_beta}, \eqn{\Sigma_\beta}{Sigma_beta} and \eqn{\lambda}{lambda} define a rather
  vague prior. The hyper-parameter \eqn{\delta \geq 2}{delta>=2} ensures the
  existence of the conditional variance. The \eqn{k}th conditional
  moment for \eqn{\epsilon_t}{e(t)} is guaranteed by setting \eqn{\delta \geq k}{delta>=k}.
  
  The Bayesian estimation of the GARCH(1,1) model with Normal
  innovations is obtained as a special case by setting \code{lambda=100}
  and \code{delta=500}. In this case, the generated values for
  \eqn{\nu}{nu} are centered around 500 which ensure approximate Normality
  for the innovations. 

  The function \code{addPriorConditions} allows to add prior conditions on the model
  parameters \eqn{\psi := (\alpha_0 \ \alpha_1 \ \beta \
    \nu)'}{psi:=(alpha0 alpha1 beta nu)'}. The
  function must return \code{TRUE} if the constraint holds and
  \code{FALSE} otherwise.

  By default, the function is:
  \preformatted{%
    addPriorConditions <- function(psi)
    {
      TRUE
    }
  }
  and therefore does not add any other constraint than the positivity of
  the parameters which are obtained through the prior distribution
  for \eqn{\psi}. 

  You simply need to modify \code{addPriorConditions} in order to add
  constraints on the model parameters \eqn{\psi}. For instance, to impose the
  covariance-stationary conditions to hold,
  i.e. \eqn{\alpha_1 + \beta < 1}, just define
  the function \code{addPriorConditions} as follows:
  \preformatted{%
    addPriorConditions <- function(psi)
    {
      psi[2] + psi[3] < 1
    }
  }
  
  Note that adding prior constraints on the model parameters
  can diminish the acceptance rate and therefore
  lead to a very inefficient sampler. This would however indicate that the
  condition is not supported by the data.

  The estimation strategy implemented in \code{bayesGARCH} is fully automatic and does not require 
  any tuning of the MCMC sampler. The generation of the Markov chains is however time
  consuming and estimating the model over several datasets on a daily basis can therefore take a significant amount
  of time. In this case, the algorithm can be easily parallelized, by running a single chain on several processors.
  Also, when the estimation is repeated over updated time series (i.e. time series with more recent 
  observations), it is wise to start the algorithm using the posterior mean or median of the parameters 
  obtained at the previous estimation step. The impact of the starting values (burn-in phase) is likely to be 
  smaller and thus the convergence faster.

  Finally, note that as any MH algorithm, the sampler can get stuck to a given value, so that the chain does not move
  anymore. However, the sampler uses Taylor-made candidate densities that are especially `constructed' at each step,
  so it is almost impossible for this MCMC sampler to get stuck at a given value for many subsequent draws.
  In the unlikely case that such ill behavior would occur, one could scale the data (to have standard deviation 1),
  or run the algorithm with different initial values or a different random seed.
  
  Please cite the package in publications. Use \code{citation("bayesGARCH")}.
}
\author{%
  David Ardia \email{david.ardia.ch@gmail.com}
}
\references{%
  Ardia, D. (2009)
  Bayesian Estimation of a Markov-Switching Threshold Asymmetric GARCH Model with Student-t Innovations.
  \emph{Econometrics Journal} \bold{12}(1), pp. 105--126. \doi{10.1111/j.1368-423X.2008.00253.x}.

  Ardia, D., Hoogerheide, L.F. (2010)
  Bayesian Estimation of the GARCH(1,1) Model with Student-t Innovations.
  \emph{The R Journal} \bold{2}(2), pp.41--47. 
  URL \url{https://journal.r-project.org/}

  Ardia, D. (2008) 
  Financial Risk Management with Bayesian Estimation of GARCH Models. 
  Lecture Notes in Economics and Mathematical Systems \bold{612}. Springer-Verlag, Berlin, Germany.
  ISBN 978-3-540-78656-6, e-ISBN 978-3-540-78657-3, doi: 10.1007/978-3-540-78657-3. 
  URL \url{http://www.springer.com/de/book/9783540786566}

  Deschamps, P.J. (2006) 
  A Flexible Prior Distribution for Markov Switching Autoregressions with Student-t Errors. 
  \emph{Journal of Econometrics} \bold{133}, pp.153--190.

  Geweke, J.F. (1993)
  Bayesian Treatment of the Independent Student-t Linear Model. 
  \emph{Journal of Applied Econometrics} \bold{8}, pp.19--40.

  Nakatsuma, T. (2000)
  Bayesian Analysis of ARMA-GARCH Models: A Markov Chain Sampling Approach. 
  \emph{Journal of Econometrics} \bold{95}(1), pp.57--69.

  Nakatsuma, T. (1998)
  A Markov-Chain Sampling Algorithm for GARCH Models'. 
  \emph{Studies in Nonlinear Dynamics and Econometrics} \bold{3}(2), pp.107--117.
}
\seealso{%
  \code{\link[fGarch:garchFit]{garchFit}} (\R package \pkg{fGarch}) for the classical
  Maximum Likelihood estimation of GARCH models.
}
\examples{%
  ## !!! INCREASE THE NUMBER OF MCMC ITERATIONS !!!

  ## LOAD DATA
  data(dem2gbp)
  y <- dem2gbp[1:750]

  ## RUN THE SAMPLER (2 chains)
  MCMC <- bayesGARCH(y, control = list(n.chain = 2, l.chain = 200))

  ## MCMC ANALYSIS (using coda)
  plot(MCMC)
  
  ## FORM THE POSTERIOR SAMPLE
  smpl <- formSmpl(MCMC, l.bi = 50)

  ## POSTERIOR STATISTICS
  summary(smpl)
  smpl <- as.matrix(smpl)
  pairs(smpl)

  ## GARCH(1,1) WITH NORMAL INNOVATIONS
  MCMC <- bayesGARCH(y, lambda = 100, delta = 500,
                     control = list(n.chain = 2, l.chain = 200))

  ## GARCH(1,1) WITH NORMAL INNOVATIONS AND 
  ## WITH COVARIANCE STATIONARITY CONDITION
  addPriorConditions <- function(psi){psi[2] + psi[3] < 1}
  MCMC <- bayesGARCH(y, lambda = 100, delta = 500,
                     control = list(n.chain = 2, l.chain = 200, 
                     addPriorConditions = addPriorConditions))
}
\keyword{models}
