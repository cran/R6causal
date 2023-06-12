#' @include R6causal.R R6causal_examples.R
NULL

#' Define structural function by a conditional probability table
#'
#' @param ycondx A data table or a data frame with the following structure
#' \itemize{
#' \item 1st column: variable to be generated, "Y"
#' \item middle columns: the parents of the the 1st column variable, "X"
#' \item last column:  the probability the case specified be the other columns, "P(Y|do(X))"
#' }
#' @param x A data table or a data frame that contains data on the variables in the
#' middle columns of \code{ycondx}, "X" and one or more columns giving data on U-variables.
#' @param Umerge_expr A character string specifying how the U-variables will
#' be combined when the value "Y" is generated, e.g. "u" or "(u1+u2)/2". The result of the
#' expression should be a random number in the interval [0,1].
#' @return A data table containing the generated  variable, "Y"
#' @examples
#' ycondx <- data.table::data.table(y =rep(c(0,1), each = 3), x=rep(1:3, 2),
#'                      prob = c(0.2,0.6,0.1,0.8,0.4,0.9))
#' x <- data.table::data.table(x = sample(1:3, 20, replace = TRUE),
#'                        uy = stats::runif(20), uy2 = stats::runif(20))
#' generate_condprob(ycondx, x, Umerge_expr = "(uy+uy2)/2")
#' @export
generate_condprob <- function(ycondx, x, Umerge_expr = NULL  ) {
  n <- nrow(x)
  # aims to avoid the possible conflict with existing names by using the prefix "R6causal_temp_".
  # initialized as NULL to avoid NOTE in R CMD check
  R6causal_temp_ID <- NULL
  R6causal_temp_prob <- NULL
  R6causal_temp_randomnum <- NULL
  R6causal_temp_lower <- NULL
  R6causal_temp_upper <- NULL
  x$R6causal_temp_ID <- 1:n
  x$R6causal_temp_dummyID <- 1
  ycondx$R6causal_temp_prob <- ycondx[,ncol(ycondx), with = FALSE]
  ycondx$R6causal_temp_dummyID <- 1
  
  if(is.null(Umerge_expr)) {
    x$R6causal_temp_randomnum <- stats::runif(n)
  } else {
    tempexpr <- paste0("x[, R6causal_temp_randomnum :=", Umerge_expr,"]")
    eval(parse(text=tempexpr))
  }
  longdt <- merge(x, ycondx, all.x = TRUE, by = intersect(names(x),names(ycondx)),
                  allow.cartesian=TRUE, sort = FALSE)
  longdt[, c("R6causal_temp_lower","R6causal_temp_upper") :=
           list(cumsum0(R6causal_temp_prob),
                cumsum(R6causal_temp_prob)),
         by = list(R6causal_temp_ID)]
  newdt <- longdt[R6causal_temp_randomnum >= R6causal_temp_lower &
                    R6causal_temp_randomnum < R6causal_temp_upper]
  resultname <- colnames(ycondx)[1]
  return(newdt[, get(resultname)])
}



#' Conduct a sequence of interventions and collect the simulated data.
#'
#' @param scm An SCM object
#' @param intervene A list where the names of the elements are the variables to be
#' intervened and the values of the elements are vectors specifying the values set
#' in the intervention
#' @param response A vector of the names of the response variables
#' @param n Size of the data to be simulated for each intervention
#' @return A list containing the values of the response variables for all intervention combinations
#' @examples
#' backdoor_experiment <- run_experiment(backdoor,
#'                                      intervene = list(x = c(0,1)),
#'                                      response = "y",
#'                                      n = 10000)
#' colMeans(backdoor_experiment$response_list$y)
#' @export
run_experiment <- function(scm, intervene, response, n) {
  response_list <- vector(mode = "list", length = length(response))
  names(response_list) <- response
  icombs <- do.call(data.table::CJ,intervene) #cross join
  resptemp <- data.table::as.data.table( matrix(NA, nrow = n, ncol = nrow(icombs)))
  for (j in 1:length(response_list)) {
    response_list[[j]] <- resptemp
  }
  for (i in 1:nrow(icombs)) {
    scmcopy <- scm$clone()
    scmcopy$intervene( names(intervene), as.numeric(icombs[i,]))
    scmcopy$simulate(n)
    for (j in 1:length(response)) {
      response_list[[ response[j] ]][,i] <- scmcopy$simdata[, get(response[j])]
    }
  }
  return(list( interventions = icombs, response_list = response_list))
}


