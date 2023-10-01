\name{CICIplus-package}
\alias{CICIplus-package}
\alias{CICIplus}
\docType{package}
\title{
Causal Inference with Continuous Interventions - Extensions for Weighted Curves
}
\description{
This package is an extension for the package \pkg{CICI}. It facilitates the estimation of counterfactual outcomes for multiple values of continuous interventions at different time points, allows plotting of causal dose-response curves and can address positivity violations through weight functions. 

It implements i) a standard g-methods approach using the \emph{sequential} \emph{g}-formula and ii) the advanced \emph{weighted} sequential \emph{g}-formula approach -- as described in Schomaker et al. (2023), see below. 

The main function of the package is currently \code{\link{sgf}}.
}

\details{
\tabular{ll}{
Package: \tab CICIplus\cr
Type: \tab Package\cr
Version: \tab 0.9.1\cr
Date: \tab 2023-08-05\cr
License: \tab GPL-2\cr
Depends: \tab R (>= 4.0), CICI\cr
Imports: \tab mgcv, SuperLearner, ggplot2, parallel, doParallel, foreach, doRNG, rngtools\cr
Suggests: \tab haldensify, hal9001 \cr
}

}
\author{
Michael Schomaker

Maintainer: Michael Schomaker <michael.schomaker@stat.uni-muenchen.de>
}

\references{
Schomaker M, McIlleron H, Denti P, Diaz I. (2023) \emph{Causal Inference for Continuous Multiple Time Point Interventions}, ArXiv e-prints: \emph{https://arxiv.org/abs/2305.06645}.
}
