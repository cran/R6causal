#' @include R6causal.R 
NULL

#' R6 Class for linear structural causal models where background variables have Gaussian distribution
#' 
#' Inherits R6 class SCM. 
#' @export
LinearGaussianSCM <- R6Class("LinearGausianSCM",
                          inherit = SCM,
                          private = list(
                            .linear_gaussian_c = NULL,
                            .linear_gaussian_B = NULL,
                            .linear_gaussian_A = NULL
                            ),
                          active = list(
                            #' @field linear_gaussian_B Matrix of structural coefficients of the observed variables
                            linear_gaussian_B = function(value) {
                              if (missing(value)) {
                                private$.linear_gaussian_B
                              } else {
                                stop("`$linear_gaussian_B` is read only", call. = FALSE)
                              }
                            },
                            #' @field linear_gaussian_A Matrix of structural coefficients of the background variables
                            linear_gaussian_A = function(value) {
                              if (missing(value)) {
                                private$.linear_gaussian_A
                              } else {
                                stop("`$linear_gaussian_A` is read only", call. = FALSE)
                              }
                            },
                            #' @field linear_gaussian_c Vector of constants in structural coefficients
                            linear_gaussian_c = function(value) {
                              if (missing(value)) {
                                private$.linear_gaussian_c
                              } else {
                                stop("`$linear_gaussian_c` is read only", call. = FALSE)
                              }
                            }),
                          public = list(
                            #' @description
                            #' Create a new linear Gaussian SCM object.
                            #' @param name Name.
                            #' @param linear_gaussian A list with the following elements:
                            #' \itemize{
                            #' \item uflist: A named list containing the functions for the background variables.
                            #' \item vnames: A vector of names of the observed variables.
                            #' \item vcoefmatrix: A matrix of coefficients for observed variables in the structural equations.
                            #' \item ucoefvector: A vector of the coefficients of dedicated error terms in the structural equations.
                            #' \item ccoefvector: A vector of constant terms in the structural equations.
                            #' \item u2coefmatrix: A matrix of the coefficients of confounding background variables in the structural equations. The number of rows equals the number of the observed variables and the number of columns equals the number of confounding background variables.
                            #' }
                            #' @param random_linear_gaussian A list with the following elements:
                            #' \itemize{
                            #' \item nv: The number of observed variables
                            #' \item edgeprob: The probability of an edge between a pair of observed variables (provide either edgeprob or avgneighbors)
                            #' \item avgneighbors: The average number of edges per a vertex (provide either edgeprob or avgneighbors)
                            #' \item u2prob: The probability of unobserved confounder between a pair of observed variables (provide either u2prob or avgu2)
                            #' \item avgu2: The average number of unobserved confounders per a vertex (provide either u2prob or avgu2)
                            #' \item vcoefdistr: A function that generates the coefficients of observed variables in the structural equations. The function must have argument 'n'.
                            #' \item ucoefdistr: A function that generates the coefficients of dedicated error terms in the structural equations. The function must have argument 'n'.
                            #' \item ccoefdistr: A function that generates the constants in the structural equations. The function must have argument 'n'.
                            #' \item u2coefdistr: A function that generates the coefficients of confounding background variables in the structural equations. The function must have argument 'n'.
                            #' }
                            #' @param rflist A named list containing the functions for missingness indicators.
                            #' @param rprefix The prefix of the missingness indicators.
                            #' @param starsuffix The suffix for variables with missing data.
                            #' @return A new `LinearGaussianSCM` object that also belongs to class `SCM`.
                            #' @examples
                            #' lgbackdoor <- LinearGaussianSCM$new("Linear Gaussian Backdoor",
                            #'                 linear_gaussian = list(
                            #'                   uflist = list(ux = function(n) {rnorm(n)},
                            #'                                 uy = function(n) {rnorm(n)},
                            #'                                 uz = function(n) {rnorm(n)}),
                            #'                   vnames = c("x","y","z"),
                            #'                   vcoefmatrix = matrix(c(0,0.4,0,0,0,0,0.6,0.8,0),3,3),
                            #'                   ucoefvector = c(1,1,1),
                            #'                   ccoefvector = c(0,0,0)))
                            #' randomlg <- LinearGaussianSCM$new("Random Linear Gaussian",
                            #'                 random_linear_gaussian = list(
                            #'                   nv = 10, 
                            #'                   edgeprob=0.5, 
                            #'                   vcoefdistr = function(n) {rnorm(n)}, 
                            #'                   ccoefdistr = function(n) {rnorm(n)}, 
                            #'                   ucoefdistr = function(n) {rnorm(n)}))
                          initialize = function(name="A linear Gaussian SCM", 
                                                linear_gaussian = NULL, 
                                                random_linear_gaussian = NULL,
                                                rflist = NULL, rprefix = "R_", starsuffix = "_md") {
                            if(is.null(linear_gaussian) && is.null(random_linear_gaussian)) {
                              stop("Either argument linear_gaussian or random_linear_gaussian is required.")
                            }
                            private$.name <- name
                            if(!is.null(linear_gaussian)) { 
                              vcoefmatrix <- linear_gaussian$vcoefmatrix
                              ucoefvector <- linear_gaussian$ucoefvector
                              ccoefvector <- linear_gaussian$ccoefvector
                              vnames <- linear_gaussian$vnames
                              nv <- length(vnames)
                              adjmatrix <- (vcoefmatrix != 0)
                              uflist <- linear_gaussian$uflist
                              unames <- names(uflist)
                              nu2 <- 0
                              if(!is.null(linear_gaussian$u2coefmatrix)) {
                                u2coefmatrix <- linear_gaussian$u2coefmatrix
                                nu2 <- ncol(u2coefmatrix)
                                u2names <- unames[(nv + 1):(nv + nu2)]
                              }
                            }
                            if(!is.null(random_linear_gaussian)) {
                              nv <- random_linear_gaussian$nv
                              edgeprob <- random_linear_gaussian$edgeprob
                              avgneighbors <- random_linear_gaussian$avgneighbors
                              if(is.null(edgeprob) & is.null(avgneighbors)) stop("random_linear_gaussian must specify either 'edgeprob' or 'avgneighbors'")
                              if(!is.null(edgeprob) & !is.null(avgneighbors)) warning("'avgneighbors' will be overrided by 'edgeprob'.")
                              if(is.null(edgeprob)) {
                                edgeprob <- min(1, 2* avgneighbors / (nv - 1)) # max number of edges is nv(nv-1)/2
                              }
                              u2prob <- random_linear_gaussian$u2prob
                              avgu2 <- random_linear_gaussian$avgu2
                              if(!is.null(u2prob) & !is.null(avgu2)) warning("'avgu2' will be overrided by 'u2prob'.")
                              if(is.null(u2prob) & !is.null(avgu2)) {
                                u2prob <- min(1, 2* avgu2 / (nv - 1)) # max number of edges is nv(nv-1)/2
                              }
                              if(is.null(u2prob)) u2prob <- 0
                              vcoefdistr <- random_linear_gaussian$vcoefdistr
                              ucoefdistr <- random_linear_gaussian$ucoefdistr
                              ccoefdistr <- random_linear_gaussian$ccoefdistr
                              u2coefdistr <- random_linear_gaussian$u2coefdistr
                              nu2 <- 0
                              if( u2prob > 0) {
                                u2matrix <- matrix( runif(nv*nv) < u2prob, nv, nv)
                                u2matrix <- u2matrix & upper.tri(u2matrix)
                                arrind <- which(u2matrix, arr.ind = TRUE)
                                nu2 <- nrow(arrind) 
                                if(nu2 > 0) {
                                  u2coefmatrix <- matrix(0, nrow = nv, ncol = nu2)
                                  for (j in 1:nu2) {
                                    u2coefmatrix[ arrind[j,1], j] <- do.call(u2coefdistr, args = list(n = 1))
                                    u2coefmatrix[ arrind[j,2], j] <- do.call(u2coefdistr, args = list(n = 1))
                                  }
                                  u2logicalmatrix <- (u2coefmatrix != 0)
                                }
                              }
                              vnames <- paste0("v",1:nv)
                              unames <- paste0("u",1:(nv+nu2))
                              if (nu2 > 0) { 
                                u2names <- unames[(nv+1):(nv+nu2)]
                              }
                              uflist <- as.list(rep("n : rnorm(n)", (nv+nu2)))
                              names(uflist) <- unames
                              adjmatrix <- matrix( runif(nv*nv) < edgeprob, nv, nv)
                              adjmatrix[lower.tri(adjmatrix)] = t(adjmatrix)[lower.tri(adjmatrix)]
                              diag(adjmatrix) <- FALSE
                              height <- runif(nv)
                              heightcomp <- outer(height,height, ">")
                              adjmatrix <-  (adjmatrix & heightcomp)
                              vcoefmatrix <- matrix( do.call(vcoefdistr, args = list(n = nv*nv)), nv, nv)
                              ucoefvector <- do.call(ucoefdistr, args = list(n = nv))
                              ccoefvector <- do.call(ccoefdistr, args = list(n = nv))
                            }
                            if(!is.null(linear_gaussian) | !is.null(random_linear_gaussian)) {
                              private$.linear_gaussian_c <- ccoefvector
                              names(private$.linear_gaussian_c) <- vnames
                              private$.linear_gaussian_B <- vcoefmatrix * adjmatrix
                              rownames(private$.linear_gaussian_B) <- vnames
                              colnames(private$.linear_gaussian_B) <- vnames
                              if (nu2 > 0) {
                                private$.linear_gaussian_A <- cbind( diag(ucoefvector), u2coefmatrix)
                              } else {
                                private$.linear_gaussian_A <- diag(ucoefvector)
                              }
                              rownames(private$.linear_gaussian_A) <- vnames
                              colnames(private$.linear_gaussian_A) <- unames
                              vflist <- vector(mode = "list", nv)
                              for (i in 1:nv) {
                                functionpartu <- paste0( ucoefvector[i], "*", unames[i])
                                if(all(!adjmatrix[i,])) {
                                  argum <- unames[i]
                                  functionpart <- paste(ccoefvector[i], "+", functionpartu)
                                } else {
                                  argum <- paste0( paste( vnames[adjmatrix[i,]], collapse = ", "),", ",unames[i])
                                  functionpartv <- paste(paste0( vcoefmatrix[i,adjmatrix[i,]], "*", vnames[adjmatrix[i,]]), collapse = "+")
                                  functionpart <- paste(ccoefvector[i], "+",  functionpartv, "+", functionpartu)
                                }
                                if( nu2 > 0 ) {
                                  if (sum(u2logicalmatrix[i,]) > 0 ) {
                                    argum <- paste0( argum, ", ", paste( u2names[u2logicalmatrix[i,]], collapse = ", "))
                                    functionpartu2 <- paste(paste0( u2coefmatrix[i,u2logicalmatrix[i,]], "*", 
                                                                    u2names[u2logicalmatrix[i,]]), collapse = " + ")
                                    functionpart <- paste(functionpart, "+", functionpartu2)
                                  }
                                }
                                vflist[[i]] <- paste(argum, ":", functionpart )
                                names(vflist) <- vnames
                              }
                              private$.is_linear_gaussian <- TRUE
                            } 
                            private$.uflist <- lapply( uflist, private$.parsefunction)
                            private$.vflist <- lapply( vflist, private$.parsefunction)
                            if(!is.null(rflist)) {
                              private$.rflist <- lapply( rflist, private$.parsefunction)
                            }
                            private$.derive_SCM()
                          }
))
                          
#' Return the mean and the covariance matrix of the conditional distribution of a linear Gaussian SCM
#'
#' @param scm A linear Gaussian SCM
#' @param situation A list with the following element:
#' \itemize{
#' \item condition : either a string that gives an SQL query ( e.g. "select x,y,z from DATA where" )
#' or a data.table consisting of the valid rows ( e.g. data.table::data.table( x = 0, y = 0))
#' }
#' @export
analytic_linear_gaussian_conditining <- function(scm, situation) {
  if(!scm$is_linear_gaussian) stop("analytic_linear_gaussian  works only for linear Gaussian SCMs (is_linear_gaussian == TRUE)")
  vnames <- scm$vnames
  unames <- scm$unames
  B <- scm$linear_gaussian_B
  A <- scm$linear_gaussian_A
  c <- scm$linear_gaussian_c
  J <- ncol(B)
  H <- ncol(A)
  I <- diag(rep(1,J))
  IB <- solve(I-B)
  muv <- IB %*% c
  rownames(muv) <- vnames
  s_vv <- IB %*% (A %*% t(A)) %*% t(IB)
  rownames(s_vv) <- vnames
  colnames(s_vv) <- vnames
  s_vu <- IB %*% A
  rownames(s_vu) <- vnames
  colnames(s_vu) <- unames
  s_uv <- t(s_vu)
  # conditioning starts
  targetnames <- names(situation$condition) 
  d <- t(situation$condition)
  freevnames <- setdiff( vnames, targetnames)
  muv1 <- muv[freevnames,,drop=FALSE]
  muv2 <- muv[targetnames,,drop=FALSE]
  s_v1v1 <- s_vv[freevnames, freevnames,drop=FALSE]
  s_v1v2 <- s_vv[freevnames, targetnames,drop=FALSE]
  s_v2v2 <- s_vv[targetnames, targetnames,drop=FALSE]
  s_v1u <- s_vu[freevnames,,drop=FALSE]
  s_v2u <- s_vu[targetnames,,drop=FALSE]
  s_block <- rbind( s_v1v2, t(s_v2u))
  inv_s_v2v2 <- solve(s_v2v2)
  mu_cond <- rbind(muv1, cbind(rep(0,H))) + s_block %*% inv_s_v2v2 %*% (d - muv2)
  rownames(mu_cond) <- c(freevnames,unames)
  s_cond <- cbind( rbind(s_v1v1, t(s_v1u)), rbind(s_v1u, diag(rep(1,H))))  - s_block %*% inv_s_v2v2 %*% t(s_block)
  colnames(s_cond) <- c(freevnames,unames)
  rownames(s_cond) <- c(freevnames,unames)
  return(list(mu = mu_cond[c(unames,freevnames),,drop=FALSE], 
              sigma = s_cond[c(unames,freevnames),c(unames,freevnames),drop=FALSE]))
}

#' Simulate data from a conditional linear Gaussian SCM
#'
#' @param scm A linear Gaussian SCM
#' @param situation A list with the following element:
#' \itemize{
#' \item condition : either a string that gives an SQL query ( e.g. "select x,y,z from DATA where" )
#' or a data.table consisting of the valid rows ( e.g. data.table::data.table( x = 0, y = 0))
#' }
#' @param n The number of rows in the data to be simulated
#' @export
analytic_linear_gaussian <- function(scm, situation, n) {
  if(!scm$is_linear_gaussian) stop("analytic_linear_gaussian  works only for linear Gaussian SCMs (is_linear_gaussian == TRUE)")
  vnames <- scm$vnames
  unames <- scm$unames
  param <- analytic_linear_gaussian_conditining(scm, situation)
  U <- MASS::mvrnorm(n = n, mu = param$mu[unames,], 
                     Sigma = param$sigma[unames,unames,drop=FALSE] )
  colnames(U) <- unames
  U <- as.data.table(U)
  simdata <- scm$simulate(n = n, fixedvars = U, return_simdata = TRUE, store_simdata = FALSE)
  return(simdata)
}
