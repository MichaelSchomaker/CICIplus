\name{sgf}
\alias{sgf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sequential g-formula for continuous multiple time point interventions
}
\description{
Estimation of counterfactual outcomes for multiple values of continuous interventions at different time points using the sequential (weighted) g-formula. 
}
\usage{
sgf(X, Anodes, Ynodes, Lnodes = NULL, Cnodes = NULL,
    abar = NULL, survivalY = FALSE, 
    SL.library = "SL.glm", SL.export = NULL,
    Yweights = NULL, calc.support = FALSE, B = 0,
    ncores = 1, verbose = TRUE, seed = NULL, prog = NULL, ...)
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
A character string of column names in \code{X} of  the time-dependent (post first treatment) variable(s).
}
  \item{Cnodes}{
A character string of column names in \code{X} of  the censoring variable(s).
}
  \item{abar}{
Numeric vector or matrix of intervention values. See Details.
}
  \item{survivalY}{
Logical. If TRUE, then Y nodes are indicators of an event, and if Y at some time point is 1, then all following should be 1.
}
  \item{SL.library}{
Either a character vector of prediction algorithms or a list containing character vectors. See details. 
}
  \item{SL.export}{
A string vector of user-written learning and screening algorithms that are not part of  \pkg{SuperLearner}, but are part of the learning library. Only required if \code{ncores>1}. See details. 
}
  \item{Yweights}{
A list of length of \code{Ynodes}, likely generated with \code{\link{calc.weights}}. 
}
  \item{calc.support}{
Logical. If \code{TRUE}, both crude and conditional support is estimated.
}

   \item{B}{
An integer specifying the number of bootstrap samples to be used, if any. 
}
  \item{ncores}{
An integer for the number of threads/cores to be used. If >1, parallelization will be utilized.
}
  \item{verbose}{
Logical. If \code{TRUE}, notes and warnings are printed.
}
  \item{seed}{
An integer specifying the seed to be used to create reproducable results for parallel computing (i.e. when ncores>1).
}
  \item{prog}{
A character specifying a path where progress should be saved (typically, when \code{ncores}>1).
}
  \item{\dots}{
Further arguments to be passed on.
}
}
\details{
The function calculates the expected counterfactual outcomes (specified under \code{Ynodes}) under the intervention \code{abar}.

If \code{abar} is a vector, then each vector component is used as the intervention value at each time point; that is, interventions which are constant over time are defined. If \code{abar} is a matrix (of size 'number interventions' x 'time points'), then each row of the length of \code{Anodes} refers to a particular time-varying intervention strategy. 

The nested iterated outcome models are fitted using super learning. The specified prediction algorithms (possibly coupled with algorithms for prior variable screening)  are passed on to package \pkg{SuperLearner}. See \code{?SuperLearner} for examples of permitted structures. Note: User-written prediction algorithms, corresponding \code{S3} prediction functions and screening algorithms need to be specified under \code{SL.export}, if parallelization is used. 

For survival settings, it is required that i) \code{survivalY=TRUE}, ii) the data are in a format where a Ynode stays 1, after it jumps to 1 and ii) after a Cnode/Ynode is 1, every variable thereafter is set to \code{NA} (except a Ynode which is already 1). See manual for an example. The package intervenes on Cnodes, i.e. calculates counterfactual outcomes under no censoring. 

If \code{calc.support=TRUE}, conditional and crude support measures (i.e., diagnostics) are calculated as described in Section 3.3.2 of Schomaker et al. (2023). 

To parallelize computations automatically, it is sufficient to set \code{ncores>1}, as appropriate. No further customization or setup is needed, everything will be done by the package. To make estimates under parallelization reproducible, use the \code{seed} argument. To watch the progress of parallelized computations, set a path in the \code{prog} argument: then, a text file reports on the progress, which is particularly useful if lengthy bootstrapping computations are required.
}
\value{
Returns an object of of \code{class} `gformula':
  \item{results}{matrix of results}
  \item{diagnostics}{list of diagnostics and weights based on the estimated support (if \code{calc.support=TRUE})}
  \item{SL.weights}{matrix of average super learner weights, at each time point}
  \item{setup}{list of chosen setup parameters}
}
\references{
Schomaker M, McIlleron H, Denti P, Diaz I. (2023) \emph{Causal Inference for Continuous Multiple Time Point Interventions}, ArXiv e-prints: \emph{https://arxiv.org/abs/2305.06645}.
}
\author{
Michael Schomaker
}

\seealso{
See \link[CICI]{gformula} for parametric g-computation and \code{\link{calc.weights}} on generating outcome weights. 
}

\examples{
data(EFV)

est <- sgf(X=EFV,
                Lnodes  = c("adherence.1","weight.1",
                            "adherence.2","weight.2",
                            "adherence.3","weight.3",
                            "adherence.4","weight.4"
                ),
                Ynodes  = c("VL.0","VL.1","VL.2","VL.3","VL.4"),
                Anodes  = c("efv.0","efv.1","efv.2","efv.3","efv.4"),
                abar=seq(0,10,1)
)
plot(est)
est

# Note: replace sgf() with gformula() for parametric g-computation
}

\keyword{sequential g-formula}
