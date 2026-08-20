# bayesGARCH 2.2.0

Update of an existing CRAN package (2.1.10 -> 2.2.0).

## Purpose of this submission

This release fixes three defects that change results. I would rather set them
out than have them noticed later.

### 1. Results depended on the units of `y`

The proposals for the scedastic parameters are built from a weighted least
squares fit whose conditioning depends on the size of `y^2` relative to
`alpha0`, whose starting value was a fixed 0.01. On a series whose scale was far
from unity, every `alpha` and `beta` move was rejected: the move rates were
exactly zero, `nu` was then dragged down to `delta`, the latent variances
collapsed, and the root solved in the `nu` update eventually left its bracket
and aborted the run. Dividing a series by 100 -- log-returns in decimals rather
than in percent -- was enough to trigger it with the default prior.

The series is now rescaled internally to unit root mean square and the draws are
mapped back. Only `alpha0` carries that scale, and the Normal prior on `alpha`
transforms exactly under the map, so the posterior explored is unchanged. The
factor is a power of two, so the rescaling and its inverse are exact in binary
floating point; a test checks that two runs differing by a power-of-two factor
agree draw for draw once the prior and the starting values are expressed in the
same units.

### 2. `addPriorConditions` did not target the constrained posterior

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

### 3. The degrees-of-freedom update could fail to bracket its root

The rejection-envelope rate was solved on the fixed interval `[1e-5, 500]`. The
root can lie above 500 -- for a sharply concentrated prior on `nu`, and formerly
through the scale defect above. The bracket now widens geometrically, and a
genuine failure to bracket reports the state rather than surfacing as `f()
values at end points not of opposite sign` from `uniroot()`.

Alongside these, inputs that cannot be represented -- prior variances whose
reciprocals overflow, a `Sigma.alpha` that cannot be inverted, a series whose
squares overflow or underflow -- are now refused with messages that say so,
instead of failing inside `solve()`, `chol()` or a foreign function call.

## Results change

Draws are no longer reproduced bit for bit from a given seed. Three things
contribute: the fixes above, the internal rescaling, and a deliberate numerical
change in which the weighted cross-products are formed from whitened regressors
through `crossprod()` and the posterior precision is inverted through its
Cholesky factor.

The posterior being targeted is unchanged. On the `dem2gbp` benchmark the
unconstrained posterior means agree with those of version 2.1.10 to within 2.0
combined Monte Carlo standard errors on every parameter. Two independent anchors
are checked in the test suite: the Normal special case reproduces the published
`dem2gbp` GARCH(1,1) benchmark, and the corrected constrained sampler is pinned
against an independent random-walk Metropolis implementation sharing no code
with the package.

There is one reverse dependency, `miscFuncs`, which suggests `bayesGARCH`. Its
`fcastGARCH()` consumes an MCMC chain matrix from `bayesGARCH()`. I checked it
against this version: the returned object, its `chain1` element, and the column
layout are unchanged, and `fcastGARCH()` produces finite forecasts as before.
Nothing downstream is affected.

## Test environments

* win-builder, R-release, R 4.6.1, x86_64-w64-mingw32, gcc 14.3.0 --
  **Status: OK** (no errors, warnings or notes; PDF and HTML manual OK)
* win-builder, R-devel, R-devel r90424 (2026-08-17), x86_64-w64-mingw32,
  gcc 14.3.0 -- **Status: OK** (no errors, warnings or notes; PDF and HTML
  manual OK)
* local: macOS 26.6.1, aarch64-apple-darwin20, R 4.5.2, Apple clang 21.0.0

## R CMD check results

No ERRORs, WARNINGs or NOTEs on win-builder, under either R-release or R-devel.

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
compile with no warnings of their own under `-Wall`, and win-builder reports
`checking whether package 'bayesGARCH' can be installed ... OK` under both
R-release and R-devel.

`checking CRAN incoming feasibility` is OK.

## Notes

* mvtnorm moved from Imports to Suggests. Its two functions used by the sampler
  were replaced by internal equivalents that perform the same operations in the
  same order; mvtnorm is now used only by tests that pin that equivalence, and
  the package no longer loads its namespace at run time.
* This version adds a testthat suite. The long-running statistical tests, which
  fit the full `dem2gbp` series, are behind `skip_on_cran()`.
