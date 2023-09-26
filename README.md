

# R/`CICIplus`


> Causal Inference for Continuous (Multiple Time Point) Interventions -- Bonus Material

**Author:** [Michael Schomaker](https://michaelschomaker.github.io/)

------------------------------------------------------------------------

## What’s `CICIplus`?

The `CICIplus` R package is an extension of the `CICI` package, which is available on CRAN. It allows estimation of counterfactual outcomes for multiple values of continuous interventions at different time points, and plotting of causal dose-response curves. This extension provides an implementation of the sequential g-formula with an integration of outcome-weights to address positivity violations, see [Schomaker, McIlleron, Denti, Diaz, 2023](https://arxiv.org/abs/2305.06645) for details. At a later stage the functions from `CICIplus` will be integrated into `CICI` . 

------------------------------------------------------------------------

## Installation

You can install this packasge from
GitHub via [`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("MichaelSchomaker/CICIplus")
```

------------------------------------------------------------------------

## Example

Examples will soon be available at [my homepage](https://michaelschomaker.github.io/project/cici/).
