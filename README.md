# bayesGARCH

Bayesian estimation of the GARCH(1,1) model with Student-t innovations

[![DOI](https://zenodo.org/badge/59887397.svg)](https://zenodo.org/badge/latestdoi/59887397)
[![CRAN](http://www.r-pkg.org/badges/version/bayesGARCH)](https://cran.r-project.org/package=bayesGARCH) 
[![Downloads](http://cranlogs.r-pkg.org/badges/bayesGARCH?color=brightgreen)](http://www.r-pkg.org/pkg/bayesGARCH)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/bayesGARCH?color=brightgreen)](http://www.r-pkg.org/pkg/bayesGARCH)

The package `bayesGARCH` ([Ardia and Hoogerheide, 2010)](https://journal.R-project.org/archive/2010-2/) implements in R 
the Bayesian estimation procedure described 
in [Ardia (2008)](http://dx.doi.org/10.1007/978-3-540-78657-3) for the GARCH(1,1) model with Student-t innovations.
The approach consists of a Metropolis-Hastings (MH) algorithm where the proposal distributions
are constructed from auxiliary ARMA processes on the squared observations. This methodology 
avoids the time-consuming and difficult task, especially for non-experts, of choosing and tuning 
a sampling algorithm. 

Please cite `bayesGARCH` in publications:

Ardia, D., Hoogerheide, L.F. (2010).  
Bayesian estimation of the GARCH(1,1) model with Student-t innovations.  
_R Journal_ **2**(2), pp.41-47.    
[https://journal.R-project.org/archive/2010-2/](https://journal.R-project.org/archive/2010-2/)    

Ardia, D. (2008).    
_Financial Risk Management with Bayesian Estimation of GARCH Models: Theory and Applications_.   
volume 612 series Lecture Notes in Economics and Mathematical Systems. Springer-Verlag, Berlin, Germany.    
[http://dx.doi.org/10.1007/978-3-540-78657-3](http://dx.doi.org/10.1007/978-3-540-78657-3)    