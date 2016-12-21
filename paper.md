---
title: 'bayesGARCH: Bayesian Estimation of the GARCH(1,1) Model with Student-t Innovations in R'
bibliography: paper.bib
date: "21 December 2016"
tags:
- GARCH
- Bayesian
- MCMC
- Student-t
- R software
authors:
- affiliation: Institute of Financial Analysis - University of Neuchâtel, Switzerland
- affiliation: Département de Finance, Assurance et Immobilier - Université Laval, Québec, Canada
  name: David Ardia
  orcid: 0000-0003-2823-782X
- affiliation: Department of Econometrics - Vrije Universiteit Amsterdam, The Netherlands
  name: Lennart F. Hoogerheide
  orcid: null
---

# Summary

The package `bayesGARCH` implements in R (@R) the Bayesian estimation procedure described 
in @Ardia2008 for the GARCH(1,1) model with Student-t innovations.
The approach consists of a Metropolis-Hastings (MH) algorithm where the proposal distributions
are constructed from auxiliary ARMA processes on the squared observations. This methodology avoids the 
time-consuming and difficult task, especially for non-experts, of choosing and tuning a sampling algorithm. 
We refer the user to @Ardia2008 and @ArdiaHoogerheide2010 for illustrations. The latest version of the 
package is available at 'https://github.com/ArdiaD/bayesGARCH'.

# References
