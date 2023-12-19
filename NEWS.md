# R6causal 0.8.2

* Fixed a bug in the interventions to missingness indicators (`rflist`).

# R6causal 0.8.1

* The content of `simdata` and `simdata_obs` revised when missing data mechanism is present. 
* Interventions can be applied also to background variables (`uflist`) and missingness indicators (`rflist`). 
* Added counterfactual identification via R package `cfid`.
* Added tikz as an output format.
* Updated the method for creating the package documentation. (I do not really understand this but it works now.)

# R6causal 0.8.0

* Added function `fairness` for fairness evaluation via counterfactual simulation.
* Added new default method "u_find" to function`counterfactual` to simulate counterfactuals for continuous and discrete variables.
* Added methods to R6 class `SCM` to add and remove vertices.
* Added methods to R6 class `SCM` to obtain parents, children, ancestors and descendants of a set of vertices.
* Added R6 class `LinearGaussianSCM` for random linear Gaussian SCMs.
* Added functions `analytic_linear_gaussian`, `analytic_linear_gaussian_conditining` to work with linear Gaussian SCMs.


# R6causal 0.7.0

* Added R6 class `ParallelWorld` to enable more general counterfactual simulations.
* Improved output control for simulations.

# R6causal 0.6.1

* A minor change due to change in the igraph package.

# R6causal 0.6.0

* The first version submitted to CRAN.
* Added a `NEWS.md` file to track changes to the package.
