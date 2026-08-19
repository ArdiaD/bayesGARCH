# bayesGARCH 2.2.0

Update of an existing CRAN package (2.1.10 -> 2.2.0).

## Purpose of this submission

This release corrects a bug that changes results for one documented use of the
package, and I would rather flag it than have it noticed later.

`bayesGARCH()` accepts a `control$addPriorConditions` function with which the
user restricts the parameter space, most commonly to the covariance-stationary
region. Two defects affected it.

First, the constraint was applied once to the whole block of parameters,
accepting or rejecting the entire sweep on its final value. Each sub-step of the
sweep is reversible with respect to the unconstrained posterior, but their
composition is only invariant and not itself reversible, so rejecting the block
does not leave the constrained posterior invariant. The chain converged to the
wrong distribution. The constraint now enters the acceptance test of each step.

Second, starting values were not checked against that constraint. Since every
step only accepts a move into the admissible region, a chain started outside it
could never enter, and returned the starting value at every iteration without
warning. Starting values are now validated and rejected with an error naming the
row and the failed condition.

Calls that did not use `addPriorConditions` were unaffected by either defect and
target the same posterior as before.

## Results change

Draws are no longer reproduced bit for bit from a given seed, for two reasons:
the bug fix above, and a deliberate numerical change in which the weighted
cross-products are formed from whitened regressors through `crossprod()` and the
posterior precision is inverted through its Cholesky factor. The posterior being
targeted is unchanged; I verified this against version 2.1.10 over long runs
(all parameters within 0.9 Monte Carlo standard errors, Kolmogorov-Smirnov
p between 0.96 and 0.99), and verified the corrected constrained sampler against
an independent random-walk Metropolis implementation sharing no code.

There is one reverse dependency, `miscFuncs`, which suggests `bayesGARCH`. Its
`fcastGARCH()` consumes an MCMC chain matrix from `bayesGARCH()`. I checked it
against this version: the returned object, its `chain1` element, and the column
layout are unchanged, and `fcastGARCH()` produces finite forecasts as before.
Nothing downstream is affected.

## Test environments

* macOS builder (mac.r-project.org), macOS arm64, R 4.6.1, clang-1700.6.3.2 --
  **Status: OK** (errors: no, warnings: no, notes: no)
* win-builder, R-devel -- **Status: OK**
* win-builder, R-release -- **Status: OK**
* local: macOS 26.6.1, aarch64-apple-darwin20, R 4.5.2, Apple clang 21.0.0

## R CMD check results

No ERRORs, WARNINGs or NOTEs on the macOS builder.

The local run reports one WARNING, which does not come from this package's
sources and does not reproduce elsewhere:

```
* checking whether package 'bayesGARCH' can be installed ... WARNING
Found the following significant warnings:
  .../R_ext/Boolean.h:62:36: warning: unknown warning group
  '-Wfixed-enum-extension', ignored [-Wunknown-warning-option]
```

It is emitted while compiling `#include <R_ext/Rdynload.h>`: R's own
`R_ext/Boolean.h` contains `#pragma clang diagnostic ignored
"-Wfixed-enum-extension"`, and the Apple clang 21 shipped with the current macOS
command line tools does not know that warning group. Both package source files
compile with no warnings of their own under `-Wall`, and the macOS builder,
which uses clang-1700, reports `checking whether package 'bayesGARCH' can be
installed ... OK`.

`checking CRAN incoming feasibility` is OK.

## Notes

* mvtnorm moved from Imports to Suggests. Its two functions used by the sampler
  were replaced by internal equivalents that perform the same operations in the
  same order; mvtnorm is now used only by tests that pin that equivalence, and
  the package no longer loads its namespace at run time.
* This version adds a testthat suite. The long-running statistical tests, which
  fit the full `dem2gbp` series, are behind `skip_on_cran()`.
