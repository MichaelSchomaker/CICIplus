\name{calc.weights}
\alias{calc.weights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate outcome weights to address positivity violations}
\description{
The weights are calculated according to formula (14) in Schomaker et al. (2023).
}
\usage{
calc.weights(X, Anodes = NULL, Ynodes = NULL, Lnodes = NULL, Cnodes = NULL,
             abar = NULL, times = length(Anodes), c = 0.01, screen = FALSE, 
             survival = FALSE, eps = 1e-10, zero = 0, 
             d.method = c("binning", "parametric", "hal_density", "hazardbinning"),
             z.method = c("density", "eps"), w.function = "gal_ga",
             for.sgf = TRUE,
             verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A data frame, following the time-ordering of the variables. 
}
  \item{Anodes}{
A character string of column names in \code{X} of the intervention variable(s).
}
  \item{Ynodes}{
A character string of column names in \code{X} of  the outcome variable(s).
}
  \item{Lnodes}{
A character string of column names in \code{X} of  \emph{all} confounders, both baseline and time-varying.
}
  \item{Cnodes}{
A character string of column names in \code{X} of  the censoring variable(s).
}
  \item{abar}{
Numeric vector or matrix of intervention values. See Details.
}
  \item{times}{
Numeric value specifying for how many time points the weights should be calculated.
}
  \item{c}{
A numeric value (or vector) specifying the threshold(s) below which the weights correspond to the density ratios, rather than 1. 
}
  \item{screen}{
Logical. If \code{TRUE}, variable screening with LASSO is performed prior to estimating the conditional densities for the weights.
}
  \item{survival}{
Logical. If \code{TRUE}, a survival setting is assumed and taken into account for model specification.
  }
  \item{eps}{
A numeric value specifying epsilon if \code{z.method="eps"}. See details. 
}
  \item{zero}{
A numeric value specifying which actual number is considered to be "zero" in the denominator. 
}
  \item{d.method}{
A string specifying which method should be used to estimate the conditional density. One of \code{"binning", "parametric", "hal_density"}.
}
  \item{z.method}{
A string specifying the method which should be used if the denominator is zero: The default is \code{"density"}, which modifies the density according to formula (14) in Schomaker et al. (2023). Alternatively, \code{"eps"} replaces the denominator with a fixed value specified under \code{eps}.
}
  \item{w.function}{
A string specifying the weight function which specifies how the conditional densities from numerator and denominator as well as c should be combined. Currently, \code{"gal_ga"} equates to formula (14) of Schomaker et al. (2023); whereas \code{"gal_ga2"} equates to the same weights, but where the cutoff for c to define positivity violations is not based on the conditional treatment density alone, but relative to the density defined in the denominator.
}
\item{for.sgf}{
 Logical. If TRUE, weights are organized such that they fit the order required for \code{\link{sgf}}.
}
  \item{verbose}{
Logical. If \code{TRUE}, notes and warnings are printed.
}
  \item{\dots}{
Further arguments to be passed on.
}
}
\details{
Calculates the outcome weights as described in formula (14) of Schomaker et al. (2023).

If \code{d.method="parametric"}, parametric conditional density estimation with generalized (additive models) is used.  Under \code{d.method="binning"} the continuous intervention is categorized into \code{length(abar)} bins, and then logistic (additive) models are used to approximate the conditional treatment density. The method \code{d.method="hal_density"} estimates the conditional treatment density non-parametrically using highly adaptive LASSO density estimation, as implemented in \pkg{haldensify}. This option is experimental so far, and may take long, especially if the sample size is large. 

In survival settings, past censoring and outcome nodes are omitted from the formulae. If censoring is present without a survival setting (e.g. Cnodes describe drop-outs and Y is a continuous outcome), then survival should be set as FALSE.
}
\value{
An object of class \code{Yweights}. This is a named list of length of \code{c}; each list entry is another list of length(number of time points); each entry is a matrix of size n times I (n=sample size; I=number of Interventions).
}
\author{
Michael Schomaker
}


\examples{
data(EFV)

w <- calc.weights(X=EFV, Lnodes  = c("sex", "metabolic", 
                                     "log_age", "NRTI" ,"weight.0",
                                     "adherence.1","weight.1",
                                     "adherence.2","weight.2",
                                     "adherence.3","weight.3",
                                     "adherence.4","weight.4"),
                        Ynodes  = c("VL.0","VL.1","VL.2","VL.3","VL.4"),
                        Anodes  = c("efv.0","efv.1","efv.2","efv.3","efv.4"),
                        d.method="binning", abar=seq(0,10,1), c=0.01)
                        
summary(w)
# w can now be used under 'Yweights' in sgf()

}

\references{
Schomaker M, McIlleron H, Denti P, Diaz I. (2023) \emph{Causal Inference for Continuous Multiple Time Point Interventions}, ArXiv e-prints: \emph{https://arxiv.org/abs/2305.06645}.
}

