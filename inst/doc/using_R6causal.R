## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(R6causal)
library(data.table)
library(stats)
data.table::setDTthreads(2)

## ----definebackdoor-----------------------------------------------------------
backdoor <- SCM$new("backdoor",
  uflist = list(
    uz = function(n) {return(runif(n))},
    ux = function(n) {return(runif(n))},
    uy = function(n) {return(runif(n))}
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


## ----definebackdoor_text------------------------------------------------------
backdoor_text <- SCM$new("backdoor",
  uflist = list(
    uz = "n : runif(n)", 
    ux = "n : runif(n)", 
    uy = "n : runif(n)" 
  ),
  vflist = list(
    z = "uz : as.numeric(uz < 0.4)", 
    x = "ux, z : as.numeric(ux < 0.2 + 0.5*z)",
    y = "uy, z, x : as.numeric(uy < 0.1 + 0.4*z + 0.4*x)" 
  )
)


## ----definebackdooralt--------------------------------------------------------
backdoor_condprob <- SCM$new("backdoor",
  uflist = list(
    uz = function(n) {return(runif(n))},
    ux = function(n) {return(runif(n))},
    uy = function(n) {return(runif(n))}
  ),
  vflist = list(
    z = function(uz) {
      return( generate_condprob( ycondx = data.table(z = c(0,1), 
                                                     prob = c(0.6,0.4)), 
                               x = data.table(uz = uz), 
                               Umerge_expr = "uz"))},
    x = function(ux, z) {
      return( generate_condprob( ycondx = data.table(x = c(0,1,0,1), 
                                                     z = c(0,0,1,1), 
                                                     prob = c(0.8,0.2,0.3,0.7)),
                                             x = data.table(z = z, ux = ux), 
                                             Umerge_expr = "ux"))},
    y = function(uy, z, x) {
      return( generate_condprob( ycondx = data.table(y= rep(c(0,1), 4), 
                                                     z = c(0,0,1,1,0,0,1,1), 
                                                     x = c(0,0,0,0,1,1,1,1), 
                                                     prob = c(0.9,0.1,0.5,0.5,
                                                              0.5,0.5,0.1,0.9)),
                                             x = data.table(z = z, x = x, uy = uy), 
                                             Umerge_expr = "uy"))}
  )
)

## ----lineargaussianbackdoor---------------------------------------------------
lgbackdoor <- LinearGaussianSCM$new("Linear Gaussian Backdoor",
                                    linear_gaussian = list(
                                      uflist = list(ux = function(n) {rnorm(n)},
                                                    uy = function(n) {rnorm(n)},
                                                    uz = function(n) {rnorm(n)}),
                                      vnames = c("x","y","z"),
                                      vcoefmatrix = matrix(c(0,0.4,0,0,0,0,0.6,0.8,0),3,3),
                                      ucoefvector = c(1,1,1),
                                      ccoefvector = c(0,0,0)))
print(lgbackdoor)

## ----randomlineargaussian-----------------------------------------------------
randomlg <- LinearGaussianSCM$new("Random Linear Gaussian",
                                  random_linear_gaussian = list(
                                  nv = 6, 
                                  edgeprob=0.5, 
                                  vcoefdistr = function(n) {rnorm(n)}, 
                                  ccoefdistr = function(n) {rnorm(n)}, 
                                  ucoefdistr = function(n) {rnorm(n)}))
print(randomlg)

## ----printbackdoor------------------------------------------------------------
backdoor

## ----plotbackdoor-------------------------------------------------------------
backdoor$plot(vertex.size = 25) # with package 'igraph'
backdoor$plot(subset = "v") # only observed variables
if (requireNamespace("qgraph", quietly = TRUE)) backdoor$plot(method = "qgraph") 
# alternative look with package 'qgraph'

## ----simulatebackdoor---------------------------------------------------------
backdoor$simulate(10)
backdoor$simdata
backdoor$simulate(8)
backdoor$simdata
backdoor_text$simulate(20)
backdoor_condprob$simulate(30)

## ----interventionbackdoor-----------------------------------------------------
backdoor_x1 <- backdoor$clone()  # making a copy
backdoor_x1$intervene("x",1) # applying the intervention
backdoor_x1$plot(method = "qgraph") # to see that arrows incoming to x are cut
backdoor_x1$simulate(10) # simulating from the intervened model
backdoor_x1$simdata

## ----interventionbackdoor2----------------------------------------------------
backdoor_yz <- backdoor$clone()  # making a copy
backdoor_yz$intervene("y", 
  function(uy, z) {return(as.numeric(uy < 0.1 + 0.8*z ))}) # making y a function of z only
backdoor_yz$plot(method = "qgraph") # to see that arrow x -> y is cut

## ----experimentbackdoor-------------------------------------------------------
backdoor_experiment <- run_experiment(backdoor, 
                                      intervene = list(x = c(0,1)), 
                                      response = "y", 
                                      n = 10000)
str(backdoor_experiment)
colMeans(backdoor_experiment$response_list$y)

## ----IDbackdoor---------------------------------------------------------------
backdoor$causal.effect(y = "y", x = "x")
backdoor$dosearch(data = "p(x,y,z)", query = "p(y|do(x))")
backdoor$cfid(gamma = cfid::conj(cfid::cf("Y",0), cfid::cf("X",0, c(Z=1))) ) 

## ----counterfactualbackdoor---------------------------------------------------
cfdata <- counterfactual(backdoor, situation = list(do = list(target = "x", ifunction = 0), 
                                                    condition = data.table( x = 0, y = 0)), 
                         target = "x", ifunction = 1, n = 100000, 
                         method = "rejection")
mean(cfdata$y)

## ----counterfactualcomparisonbackdoor-----------------------------------------
backdoor_x1$simulate(100000)
mean(backdoor_x1$simdata$y)

## ----parallelworldbackdoor----------------------------------------------------
backdoor_parallel <- ParallelWorld$new(
                         backdoor,
                        dolist=list(
                             list(target = "x", 
                                 ifunction = 0),
                             list(target = list("z","x"), 
                                  ifunction = list(1,0))
                         )
 )
backdoor_parallel
if (requireNamespace("qgraph", quietly = TRUE)) backdoor_parallel$plot(method = "qgraph") 

## ----parallelworldcounterfactual----------------------------------------------
cfdata <- counterfactual(backdoor_parallel,
                         situation = list(
                            do = NULL,
                            condition = data.table::data.table( y = 0, y_1 = 0, y_2 = 0)),
                         target = "x",
                         ifunction = 1,
                         n = 100000,
                         method = "rejection")
mean(cfdata$y)

## ----parallelworldbackdoor2---------------------------------------------------
backdoor_parallel2 <- ParallelWorld$new(
                         backdoor,
                        dolist=list(
                             list(target = "x", 
                                 ifunction = 0),
                             list(target = list("z","x"), 
                                  ifunction = list(1,0)),
                              list(target = "x", 
                                 ifunction = 1)
                         )
 )
cfdata <- counterfactual(backdoor_parallel2,
                         situation = list(
                            do = NULL,
                            condition = data.table::data.table( y = 0, y_1 = 0, y_2 = 0)),
                         n = 100000, 
                         method = "rejection")
mean(cfdata$y_3)

## ----definebackdoor_md--------------------------------------------------------
backdoor_md <- SCM$new("backdoor_md",
                       uflist = list(
                         uz = "n : runif(n)",
                         ux = "n : runif(n)",
                         uy = "n : runif(n)",
                         urz = "n : runif(n)",
                         urx = "n : runif(n)",
                         ury = "n : runif(n)"
                       ),
                       vflist = list(
                         z = "uz : as.numeric(uz < 0.4)",
                         x = "ux, z : as.numeric(ux < 0.2 + 0.5*z)",
                         y = "uy, z, x : as.numeric(uy < 0.1 + 0.4*z + 0.4*x)"
                       ),
                       rflist = list(
                         z = "urz : as.numeric( urz < 0.9)",
                         x = "urx, z : as.numeric( (urx + z)/2 < 0.9)",
                         y = "ury, z : as.numeric( (ury + z)/2 < 0.9)"
                       ),
                       rprefix = "r_"
)

## ----plotbackdoor_md----------------------------------------------------------
backdoor_md$plot(vertex.size = 25, edge.arrow.size=0.5) # with package 'igraph'
backdoor_md$plot(subset = "v") # only observed variables a
if (!requireNamespace("qgraph", quietly = TRUE)) backdoor_md$plot(method = "qgraph") 
# alternative look with package 'qgraph'

## ----backdoor_md_simulate-----------------------------------------------------
backdoor_md$simulate(100)
summary(backdoor_md$simdata)
summary(backdoor_md$simdata_obs)

## ----backdoor_md_simulate2----------------------------------------------------
backdoor_md$simulate(100, fixedvars = c("x","y","z","ux","uy","uz"))
summary(backdoor_md$simdata)
summary(backdoor_md$simdata_obs)

## ----backdoor_md_dosearch-----------------------------------------------------
backdoor_md$dosearch(data = "p(x*,y*,z*,r_x,r_y,r_z)", query = "p(y|do(x))")

