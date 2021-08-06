
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R6causal

<!-- badges: start -->
<!-- badges: end -->

The goal of R6causal is to make it easy to work with structural causal
models. The R6 class ‘SCM’ contains methods for 1) defining a structural
causal model via functions, text or conditional probability tables, 2)
printing basic information on the model, 3) plotting the graph for the
model using packages ‘igraph’ or ‘qgraph’, 4) simulating data from the
model, 5) applying an intervention, 6) checking the identifiability of a
query using the R packages ‘causaleffect’ and ‘dosearch’, 7) defining
the missing data mechanism, 8) simulating incomplete data from the model
according to the specified missing data mechanism and 9) checking the
identifiability in a missing data problem using the R package
‘dosearch’. In addition, there are functions for running experiments and
doing counterfactual inference using simulation.

## Installation

You can install the released version of R6causal from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("R6causal")
```

## Example

This is a basic example which shows you how define a structural causal
model, plot the graph and simulate data from the model.

``` r
library(R6causal)
backdoor <- SCM$new("backdoor",
  uflist = list(
    uz = function(n) {return(stats::runif(n))},
    ux = function(n) {return(stats::runif(n))},
    uy = function(n) {return(stats::runif(n))}
  ),
  vflist = list(
    z = function(uz) {
      return(as.numeric(uz < 0.4))},
    x = function(ux, z) {
      return(as.numeric(ux < 0.2 + 0.5*z))},
    y = function(uy, z, x) {
      return(as.numeric(uy < 0.1 + 0.4*z + 0.4*x))}
  )
)
#backdoor$plot(vertex.size = 25)
backdoor$simulate(10)
backdoor$simdata
#>             uz         ux          uy z x y
#>  1: 0.95831862 0.39196794 0.321062576 0 0 0
#>  2: 0.43000353 0.02233175 0.617388078 0 1 0
#>  3: 0.48561035 0.91441259 0.001838702 0 0 1
#>  4: 0.06472618 0.82086014 0.390221544 1 0 1
#>  5: 0.36554423 0.51432070 0.203865336 1 1 1
#>  6: 0.75471971 0.57672430 0.391862174 0 0 0
#>  7: 0.75034180 0.85961544 0.652490875 0 0 0
#>  8: 0.46957831 0.20078767 0.608333027 0 0 0
#>  9: 0.09367406 0.29394365 0.488680623 1 1 1
#> 10: 0.92747329 0.84545162 0.230997331 0 0 0
```
