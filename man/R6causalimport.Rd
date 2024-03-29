% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R6causal.R
\name{R6causalimport}
\alias{R6causalimport}
\title{R6causal: R6 class for structural causal models}
\description{
Package R6causal implements an R6 class for structural causal models (SCM)
with latent variables and missing data mechanism.  The class contains methods for
1) defining a structural causal model via functions, text or conditional probability tables,
2) printing basic information on the model,
3) plotting the graph for the model using packages `igraph` or `qgraph`,
4) simulating data from the model, 5) applying an intervention,
6) checking the identifiability of a query using the R packages `causaleffect` and `dosearch`,
7) defining the missing data mechanism,
8) simulating incomplete data from the model according to the specified missing data mechanism and
9) checking the identifiability in a missing data problem using the R package `dosearch`.
In addition, there are functions for running experiments and doing counterfactual inference using simulation.
}
\references{
J. Pearl (2009). Causality, 2nd edition, Cambridge University Press.
}
