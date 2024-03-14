#' @include R6causal.R R6causal_examples.R
NULL

#' R6 Class for parallel world models
#' 
#' Inherits R6 class SCM. 
#' @export
ParallelWorld <- R6Class("ParallelWorld",
                         inherit = SCM,
                         private = list(
                           .num_worlds = NULL,
                           .worldnames = NULL,
                           .worldsuffix = NULL,
                           .originalscm = NULL,
                           .dolist = NULL
                         ),
                         active = list(
                           #' @field num_worlds Number of parallel worlds.
                           num_worlds = function(value) {
                             if (missing(value)) {
                               private$.num_worlds
                             } else {
                               stop("`$num_worlds` is read only", call. = FALSE)
                             }
                           },
                           #' @field worldnames Names of parallel worlds.
                           worldnames = function(value) {
                             if (missing(value)) {
                               private$.worldnames
                             } else {
                               private$.worldnames <- value
                             }
                           },
                           #' @field worldsuffix Suffix used for parallel world variables.
                           worldsuffix = function(value) {
                             if (missing(value)) {
                               private$.worldsuffix
                             } else {
                               stop("`$worldsuffix` is read only", call. = FALSE)
                             }
                           },
                           #' @field originalscm SCM from which the parallel worlds are derived.
                           originalscm = function(value) {
                             if (missing(value)) {
                               private$.originalscm
                             } else {
                               stop("`$originalscm` is read only", call. = FALSE)
                             }
                           },
                           #' @field dolist List containing the interventions for each world.
                           dolist = function(value) {
                             if (missing(value)) {
                               private$.dolist
                             } else {
                               stop("`$dolist` is read only", call. = FALSE)
                             }
                           }
                         ),
                         public = list(
                           #' @description
                           #' Create a new ParallelWorld object from an SCM object.
                           #' @param scm An SCM object.
                           #' @param dolist A list containing the interventions for each world. Each element 
                           #' of the list has the fields:
                           #' \itemize{
                           #' \item target: a vector of variable names that specify the target
                           #' variable(s) of the counterfactual intervention.
                           #' \item ifunction: a list of functions for the counterfactual intervention.
                           #' }
                           #' @param worldnames A character vector giving the names of the parallel worlds.
                           #' @param worldsuffix A text giving the suffix used for parallel world variables 
                           #' before the world number. Defaults to "_" and the worlds have then suffixes 
                           #' "_1", "_2", "_3", ...
                           #' @return A new `ParallelWorld` object that also belongs to class `SCM`.
                           #' @examples
                           #' backdoor_parallel <- ParallelWorld$new(
                           #'                         backdoor,
                           #'                         dolist=list(
                           #'                             list(target = "x", 
                           #'                                  ifunction = 0),
                           #'                             list(target = list("z","x"), 
                           #'                                  ifunction = list(1,0))
                           #'                         )
                           #' )
                           #' backdoor_parallel 
                           #' plot(backdoor_parallel)
                           initialize = function(scm, dolist, worldnames = NULL, 
                                                 worldsuffix = "_") {
                             private$.originalscm <- scm$clone(deep = TRUE) #.originalscm must be read only
                             private$.dolist <- dolist
                             private$.name <- scm$name
                             if(!is.null(worldnames)) {
                               private$.worldnames <- worldnames 
                             } else {
                               private$.worldnames <- 1:length(dolist)
                             }
                             private$.worldsuffix <- worldsuffix
                             private$.uflist <-  scm$uflist
                             private$.rflist <-  scm$rflist
                             private$.rprefix <- scm$rprefix
                             private$.starsuffix <- scm$starsuffix
                             vflist <- scm$vflist
                             for(j in 1:length(dolist)) {
                               twin <- scm$clone()
                               twin$intervene(dolist[[j]]$target, dolist[[j]]$ifunction) 
                               newnames <- paste0(names(twin$vflist), worldsuffix, j)
                               names(newnames) <- names(twin$vflist)
                               for(k in 1:length(twin$vflist)) {
                                 newfunction <- rename_arguments(twin$vflist[[k]], newnames)
                                 vflist <- c(vflist, newfunction) 
                                 names(vflist)[length(vflist)] <- newnames[k] 
                               }
                             }
                             private$.vflist <- lapply( vflist, private$.parsefunction)
                             
                             private$.num_u <- length(private$.uflist)
                             private$.num_v <- length(private$.vflist)
                             private$.num_r <- length(private$.rflist)
                             private$.num_uv <- private$.num_u + private$.num_v
                             private$.num_vr <- private$.num_v + private$.num_r
                             private$.num_uvr <- private$.num_u + private$.num_v + private$.num_r
                             private$.vnames <- names(private$.vflist)
                             private$.unames <- names(private$.uflist)
                             private$.uvnames <- c( private$.unames, private$.vnames)
                             private$.derive_graph()
                             if( !is.null(private$.rflist)) {
                               private$.rprefix <- private$.rprefix
                               private$.starsuffix <- private$.starsuffix
                               private$.rnames_target  <- names(private$.rflist)
                               private$.rnames_prefix  <- paste0( private$.rprefix, 
                                                                  names(private$.rflist))
                               private$.rnames <- private$.rnames_prefix
                               private$.rflist_prefix <- private$.rflist
                               names(private$.rflist_prefix) <- private$.rnames
                               private$.rmapping_from_prefix <- as.list(names(private$.rflist))
                               names(private$.rmapping_from_prefix) <- private$.rnames
                               private$.rmapping_to_prefix <- as.list(private$.rnames)
                               names(private$.rmapping_to_prefix) <- names(private$.rflist)
                               private$.vrflist <- c(private$.vflist, private$.rflist_prefix)
                               private$.vrnames <- names(private$.vrflist)
                               private$.uvrnames <- c(private$.unames, private$.vrnames)
                               private$.derive_graph_md()
                             }
                           }
                         )
)


u_find_onecondition <- function(twin, situation, fixedvars = NULL, n, control, 
                                weightfunction = stats::dnorm, minu = -10, maxu = 10) {
  minu_ini <- minu
  maxu_ini <- maxu
  max_iterations <- control$max_iterations
  batchsize <- control$batchsize 
  if(!is.null(fixedvars)) batchsize <- min(batchsize, nrow(twin$simdata))
  sampling_replace <- control$sampling_replace
  diffmaxuminu <- maxu_ini - minu_ini
  batchi <- 1
  targetname <- names(situation$condition) #One condition only
  targetvalue <- as.numeric(situation$condition)
  targetfunction <- twin$vflist[[targetname]]
  targetargs <- names(formals(targetfunction))
  targetu <- twin$dedicated_u[[targetname]]
  minu <- rep(minu_ini,batchsize)
  maxu <- rep(maxu_ini,batchsize)
  twin$simulate(n = batchsize, fixedvars = fixedvars)
  twin$simdata[,c(targetu) := maxu] 
  maxuvalue <- do.call("targetfunction", twin$simdata[, targetargs, with = FALSE])
  twin$simdata[,c(targetu) := minu] 
  minuvalue <- do.call("targetfunction", twin$simdata[, targetargs, with = FALSE])
  signmaxmin <- sign( sign(targetvalue - maxuvalue) - sign(targetvalue - minuvalue))
  notok_ind <- (signmaxmin == 0)
  if(sum(!notok_ind) == 0) {
    stop(paste("Could not find any valid value of dedicated error term when 
               processing condition",targetname,"=",targetvalue))
  }
  signmaxmin <- signmaxmin[!notok_ind]
  twin$simdata[,c(targetu) := 0] 
  for (j in seq(-2, -max_iterations, by = -1)) {
    currentvalue <- do.call("targetfunction", twin$simdata[!notok_ind, targetargs, with = FALSE])
    twin$simdata[!notok_ind, c(targetu) := twin$simdata[!notok_ind, targetu, with = FALSE] - 
                   2^j * sign( targetvalue - currentvalue) * 
                   signmaxmin * diffmaxuminu]
  }
  uprobs <- rep(0,n)
  #uprobs[!notok_ind] <- stats::dnorm(unlist(twin$simdata[!notok_ind, targetu, with = FALSE], use.names = FALSE))
  args <- twin$simdata[!notok_ind, targetu, with = FALSE]
  names(args) <- names(formals(weightfunction)[1])
  uprobs[!notok_ind] <- do.call("weightfunction", args)
  valid_ind <- sample.int(length(uprobs), size = n, 
                          replace = sampling_replace, prob = uprobs)
  twin$simdata <- twin$simdata[valid_ind,]
  twin$simulate(n = n, fixedvars = twin$unames)
  return(NULL)
}

makelist <- function(x,targetnames) {
  if(length(x) == 1) {
    xlist <- vector(length(targetnames), mode = "list")
    names(xlist) <- targetnames
    for(k in 1:length(targetnames)) {
      xlist[[k]] <- x
    }
  } else {
    xlist <- x 
  }
  return(xlist)
}

u_find <- function(twin, situation, n, control) {
  batchsize <- control$batchsize 
  nonunique_jittersd <- control$nonunique_jittersd
  condition_type <- situation$condition_type
  weightfunction <- control$weightfunction
  minu <- control$minu
  maxu <- control$maxu
  targetnames <- names(situation$condition) 
  targetnames <- intersect(twin$toporderv, targetnames) #Sorts the names in topological order (order not specified in the documentation of intersect)
  targetsituation <- vector(length(targetnames), mode = "list")
  targetancestors <- vector(length(targetnames), mode = "list")
  for(k in 1:length(targetnames)) {
    targetsituation[[k]] <- list(condition = subset(situation$condition, select = targetnames[k] )) 
    targetancestors[[k]] <- names(unlist(igraph::ego(graph = twin$igraph,
                                                     order = length(twin$toporder),
                                                     nodes = targetnames[k] ,mode="in")))
  }
  variables_processed <- NULL
  if(is.null(condition_type)) {
    stop(" situation$condition_type must be specified as a character vector giving the type ('continuous' or 'discrete') of every variable situation$condition." )
  }
  weightfunctionlist <- makelist(weightfunction, targetnames)
  minulist <- makelist(minu, targetnames)
  maxulist <- makelist(maxu, targetnames)
  for(k in 1:length(targetnames)) {
    #print(targetnames[[k]])
    matched_condition_type <- match.arg(condition_type[[ targetnames[[k]] ]],
                                        c("continuous","discrete","categorical"))    
    if(k == 1) {
      if( identical(matched_condition_type, "continuous")) {
        u_find_onecondition(twin, situation = targetsituation[[k]], 
                            fixedvars = NULL, n = batchsize, control,
                            weightfunction = weightfunctionlist[[ targetnames[[k]] ]],
                            minu = minulist[[ targetnames[[k]]]], 
                            maxu = maxulist[[ targetnames[[k]]]])
      } 
      if( matched_condition_type %in% c("discrete","categorical")) {
        discrete_sampling_onecondition(twin, situation = targetsituation[[k]], 
                                        fixedvars = NULL, n = batchsize, control)
      }
    } else {
      if( identical(matched_condition_type, "continuous") ) {
        u_find_onecondition(twin, situation = targetsituation[[k]], 
                            fixedvars = union(variables_processed, twin$de(targetnames[[k]])), 
                            n = nrow(twin$simdata), control,
                            weightfunction = weightfunctionlist[[ targetnames[[k]] ]],
                            minu = minulist[[ targetnames[[k]]]], 
                            maxu = maxulist[[ targetnames[[k]]]])
      }
      if( matched_condition_type %in% c("discrete","categorical")) {
        discrete_sampling_onecondition(twin, situation = targetsituation[[k]], 
                                        fixedvars = union(variables_processed, twin$de(targetnames[[k]])), 
                                        n = nrow(twin$simdata), control)
      }
    }
    variables_processed <- unique(c(variables_processed, targetancestors[[k]]))
  }
  if(nrow(twin$simdata) < n) {
    warning(paste("The number of rows in the simulated data is only",nrow(twin$simdata)))
  }
  cf_data <- twin$simdata[1:min(n,nrow(twin$simdata)),]
  if(!is.null(nonunique_jittersd)) {
    targetancestors_all <- unique(unlist(targetancestors))
    duplicated_ind <- duplicated(cf_data[, targetancestors_all, with = FALSE])
    nduplicated <- sum(duplicated_ind)
    if(nduplicated > 0) {
      cf_data[duplicated_ind, targetancestors_all] <- 
        cf_data[duplicated_ind, targetancestors_all, with = FALSE] + 
        matrix( stats::rnorm( nduplicated * length(targetancestors_all), 0, nonunique_jittersd), 
                nrow = nduplicated, ncol = length(targetancestors_all))
      cf_data <- twin$simulate(n = nrow(cf_data), 
                               fixedvars = subset(cf_data, select = twin$unames),
                               store_simdata = FALSE, return_simdata = TRUE)
    }
  }
  return(cf_data)
}


# discrete_sampling_onecondition <- function(twin, situation, fixedvars = NULL, n, control) {
#   twin$simulate(n = n, fixedvars = fixedvars)
#   validdata <- merge( twin$simdata, situation$condition, by = names(situation$condition), allow.cartesian = TRUE)
#   if(nrow(validdata) == 0) {
#     stop(paste("Could not find any valid values of background variables when processing condition",names(situation$condition),"=",situation$condition))
#   }
#   valid_ind <- sample.int(nrow(validdata), size = n, replace = TRUE)
#   twin$simdata <- validdata[valid_ind,]
#   twin$simulate(n = n, fixedvars = twin$unames)
#   return(twin$simdata)
# }

discrete_sampling_onecondition <- function(twin, situation, fixedvars = NULL, n, control) {
  # \itemize{
  # \item batchsize A scalar, the number of rows to be simulated for a batch (rejection sampling)
  # \item maxbatchs A scalar, the maximum number of batchs (rejection sampling)
  # }
  batchsize <- control$batchsize
  maxbatchs <- control$maxbatchs
  if(is.null(batchsize )) batchsize <- n
  if(is.null(maxbatchs )) maxbatchs <- 1
  twin$simulate(n = batchsize, fixedvars = fixedvars)
  cf_data <- merge( twin$simdata, situation$condition, by = names(situation$condition), allow.cartesian = TRUE)
  batchi <- 2
  while( nrow(cf_data) < n & batchi <= maxbatchs) {
    twin$simulate(n = batchsize, fixedvars = fixedvars)
    validdata <- merge( twin$simdata, situation$condition, by = names(situation$condition), allow.cartesian = TRUE)
    cf_data <- rbind(cf_data, validdata)
    batchi <- batchi + 1
  }
  if(nrow(cf_data) == 0) {
    stop(paste("Could not find any valid values of background variables when processing condition",names(situation$condition),"=",situation$condition))
  }
  valid_ind <- sample.int(nrow(cf_data), size = n, replace = TRUE)
  twin$simdata <- cf_data[valid_ind,]
  twin$simulate(n = n, fixedvars = twin$unames)
  return(twin$simdata)
  #return(cf_data)
}


# rejection_sampling_onecondition <- function(twin, situation, fixedvars = NULL, n, control) {
#   # \itemize{
#   # \item batchsize A scalar, the number of rows to be simulated for a batch (rejection sampling)
#   # \item maxbatchs A scalar, the maximum number of batchs (rejection sampling)
#   # }
#   batchsize <- control$batchsize 
#   maxbatchs <- control$maxbatchs
#   twin$simulate(1)
#   cf_data <- twin$simdata[0]
#   batchi <- 1
#   while( nrow(cf_data) < n & batchi <= maxbatchs) {
#     twin$simulate(n = batchsize, fixedvars = fixedvars)
#     if( identical( class(situation$condition), "character")) {
#       if (!requireNamespace("sqldf", quietly = TRUE)) {
#         stop("Package \"sqldf\" needed when the condition is given as SQL. Please install it.",
#              call. = FALSE)
#       }
#       DATA <- twin$simdata
#       validdata <- sqldf::sqldf(situation$condition)
#     } else {
#       validdata <- merge( twin$simdata, situation$condition, by = names(situation$condition), allow.cartesian = TRUE)
#     }
#     cf_data <- rbind(cf_data, validdata)
#     batchi <- batchi + 1
#   }
#   if(nrow(cf_data) == 0) {
#     stop(paste("Could not find any valid values of background variables when processing condition",names(situation$condition),"=",situation$condition))
#   }
#   return(cf_data)
# }


scm_rejection_sampling <- function(twin, situation, n, control) {
  # \itemize{
  # \item batchsize A scalar, the number of rows to be simulated for a batch (rejection sampling)
  # \item maxbatchs A scalar, the maximum number of batchs (rejection sampling)
  # }
  batchsize <- control$batchsize 
  maxbatchs <- control$maxbatchs
  twin$simulate(1)
  cf_data <- twin$simdata[0]
  batchi <- 1
  while( nrow(cf_data) < n & batchi <= maxbatchs) {
    twin$simulate(n = batchsize)
    if( identical( class(situation$condition), "character")) {
      if (!requireNamespace("sqldf", quietly = TRUE)) {
        stop("Package \"sqldf\" needed when the condition is given as SQL. Please install it.",
             call. = FALSE)
      }
      DATA <- twin$simdata
      validdata <- sqldf::sqldf(situation$condition)
    } else {
      validdata <- merge( twin$simdata, situation$condition, by = names(situation$condition), allow.cartesian = TRUE)
    }
    cf_data <- rbind(cf_data, validdata)
    batchi <- batchi + 1
  }
  return(cf_data)
}


#' Counterfactual inference via simulation
#'
#' @param scm An SCM object
#' @param situation A list or a character string. The list has the following elements:
#' \itemize{
#' \item do : NULL or a list containing named elements 'target' and 'ifunction' that
#' specify the intervention carried out in the situation
#' \item dolist : NULL or a list of lists containing named elements 'target' and 
#' 'ifunction' that specify the intervention carried out in each parallel world
#' \item condition : either a string that gives an SQL query ( e.g. "select x,y,z from DATA where" )
#' or a data.table consisting of the valid rows ( e.g. data.table::data.table( x = 0, y = 0))
#' \item condition_type : (required only if method == "u_find") A character vector giving the type ("continuous" or "discrete") of every variable in \code{situation$condition}
#' }
#' @param n The number of rows in the data to be simulated
#' @param method The simulation method, "u_find", "rejection" or "analytic_linear_gaussian"
#' @param target NULL or a vector of variable names that specify the target
#' variable(s) of the counterfactual intervention.
#' @param ifunction NULL or a list of functions for the counterfactual intervention.
#' @param returnscm A logical, should the internally created twin SCM or parallel 
#' world SCM returned?
#' @param control List of parameters to be passed to the simulation method:
#' \itemize{
#' \item batchsize: (u_find, rejection) The size of data from n observations are resampled (default n)
#' \item max_iterations: (u_find) The maximum number of iterations for the binary search (default 50)
#' \item minu: (u_find) A scalar or a named list that specifies the lower starting value for the binary search (default -10)
#' \item maxu: (u_find) A scalar or a named list that specifies the upper starting value for the binary search (default 10)
#' \item sampling_replace: (u_find) Logical, resampling with replacement? (default TRUE)
#' \item nonunique_jittersd: (u_find) Standard deviation of the noise to be added to the output (default NULL meaning no noise)
#' \item maxbatchs: (u_find, rejection) The maximum number of batches for rejection sampling (for discrete variables)
#' \item weightfunction: (u_find) A function or a named list of functions to be applied to dedicated error terms to obtain the resampling weights (default stats::dnorm)
#' }
#' @return A data table representing the situation after the counterfactual intervention
#' @examples
#' cfdata <- counterfactual(backdoor,
#'                          situation = list(
#'                              do = list(target = "x", ifunction = 0),
#'                              condition = data.table::data.table( x = 0, y = 0)),
#'                          target = "x",
#'                          ifunction = 1,
#'                          method = "rejection",
#'                          n = 1000)
#' mean(cfdata$y)
#' 
#' backdoor_parallel <- ParallelWorld$new(backdoor,
#'                                        dolist=list(
#'                                          list(target = "x", ifunction = 0),
#'                                          list(target = list("z","x"), ifunction = list(1,0))
#'                                        )
#' )
#' cfdata2 <- counterfactual(backdoor_parallel,
#'                          situation = list(
#'                              do = NULL,
#'                              condition = data.table::data.table( y = 0, y_1 = 0, y_2 = 0)),
#'                          target = "x",
#'                          ifunction = 1,
#'                          method = "rejection",
#'                          n = 1000)
#' mean(cfdata2$y)
#' @export
counterfactual <- function(scm, situation, n, target = NULL, ifunction = NULL, method = NULL, returnscm = FALSE, control = NULL) {
  control_defaults <- list(batchsize = n, 
                           max_iterations = 50,
                           minu = -10,
                           maxu = 10,
                           sampling_replace = TRUE,
                           nonunique_jittersd = NULL,
                           maxbatchs = 100,
                           weightfunction = stats::dnorm
  )
  control_new <- control_defaults
  if(is.null(control)) {
    control <- control_defaults   
  } else {
    for(j in 1:length(control)) {
      ind <- try(match.arg( names(control)[[j]], names(control_defaults)))
      if(!inherits(ind, "try-error")) control_new[[ind]] <- control[[ind]]
    }
    control <- control_new
  }
  if( !(method %in% c("rejection","analytic_linear_gaussian","u_find"))) {
    stop(" 'method' must have one of the values in c('u_find', 'rejection','analytic_linear_gaussian')")
  }  
  # if( !(method %in% c("rejection","analytic_linear_gaussian","u_find") || 
  #       control$method %in% c("rejection","analytic_linear_gaussian","u_find"))) {
  #  stop("Either 'method' or 'control$method' must have one of the values in c('rejection','analytic_linear_gaussian','u_find')")
  # }
  # if(is.null(control$method)) {
  #   control$method <- method
  # }
  # if( !(identical(method,control$method))) {
  #   stop(" 'method' and 'control$method' differ from each other. Please specify only one of these alternatives." )
  # }
  if(!is.null(situation$dolist)) {
    twin <- ParallelWorld$new(scm, situation$dolist)     
  } else {
    twin <- scm$clone()
    if(!is.null( situation$do)) {
      twin$intervene(situation$do$target, situation$do$ifunction)
    }
  }
  matched_method<- match.arg( method, c("rejection","analytic_linear_gaussian","u_find"))  
  if( identical(matched_method, "rejection")) {
    cf_data <- scm_rejection_sampling(twin, situation, n, control)
  } 
  if( identical(matched_method, "analytic_linear_gaussian")) {
    cf_data <- analytic_linear_gaussian(twin, situation, n)
  } 
  if( identical(matched_method, "u_find")) {
    cf_data <- u_find(twin, situation, n, control)
  } 
  truen <- min(n, nrow(cf_data))
  if(!identical(nrow(twin$simdata), truen)) twin$simulate(truen)
  twin$simdata <- cf_data[1:truen, ]
  if(!is.null(target)) {
    twin$intervene(target, ifunction)
    fixedvars <- names( twin$simdata)
    fixedvars <- setdiff( fixedvars, descendants(target, scm$igraph) )
    twin$simulate(truen, fixedvars =  fixedvars)
  } 
  if(returnscm) {
    return(twin) 
  } else {
    return( twin$simdata)
  }
}

#' Checking fairness of a prediction via counterfactual simulation
#'
#' @param modellist A list of model objects that have a predict method or a list of functions that return predictions
#' @param scm An SCM object
#' @param sensitive A character vector of the names of sensitive variables
#' @param condition A data.table consisting of the valid rows ( e.g. data.table::data.table( x = 0, y = 0))
#' @param condition_type (required only if method == "u_find") A character vector giving the type ("continuous" or "discrete") of every variable in \code{condition}
#' @param parents A character vector of the names of variables that remain fixed
#' @param n The number of rows in the data to be simulated by \code{counterfactual}
#' @param sens_values A data.table specifying the combinations of the values of sensitive variables to be considered (default NULL meaning the all possible combinations of the values of sensitive variables)
#' @param modeltype "predict" (default) or "function" depending on the type \code{modellist}
#' @param method The simulation method, "u_find", "rejection" or "analytic_linear_gaussian"
#' @param control List of parameters to be passed to the simulation method, see \code{counterfactual}.
#' @param ... Other arguments passed to \code{predict} or to the prediction functions.
#' @return A list containing a data table for element of \code{modellist}. Each data table contains the predicted values after counterfactual interventions on the sensitive variables.
#' @examples
#' trainingd <- backdoor$simulate(10000, return_simdata = TRUE)
#' newd <- backdoor$simulate(100, return_simdata = TRUE)
#' vnames <- backdoor$vnames
#' m1 <- lm(y ~  x + z, data = trainingd)
#' m2 <- lm(y ~  z, data = trainingd)
#' fairlist <- fairness(modellist = list(m1,m2),
#'                      scm = backdoor,
#'                      sensitive = c("x"),
#'                      sens_values = data.table::data.table(x=c(0,1)),
#'                      condition = newd[1,c("x","y")],
#'                      condition_type = list(x = "cont",
#'                                            z = "cont",
#'                                            y = "cont"),
#'                      parents = NULL,
#'                      n = 20,
#'                      modeltype = "predict",
#'                      method = "u_find")
#' @export
fairness <- function(modellist, scm, sensitive, condition, condition_type, 
                     parents, n, sens_values = NULL, modeltype = "predict",
                     method,  control = NULL, ...) {
  if(length(modeltype) == 1) modeltype <- rep(modeltype, length(modellist))
  cfscm <- counterfactual(scm, situation = list(condition = condition, condition_type = condition_type), 
                          n = n, method = method, returnscm = TRUE, control = control)
  if(is.null(sens_values)) {
    sens_values <- unique(cfscm$simdata[, sensitive, with = FALSE])
  }
  sens_values_text <- apply(sens_values,1,paste0,collapse="_")
  predictlist <- vector(length(modellist), mode = "list")
  for (j in 1:length(modellist)) {
    predictlist[[j]] <- as.data.table(matrix(NA, nrow = n, ncol = nrow(sens_values)))
    colnames(predictlist[[j]]) <- sens_values_text
  }
  for (i in 1:nrow(sens_values)) {
    twin <- cfscm$clone()
    twin$intervene(target = sensitive, ifunction = sens_values[i,])
    if(!is.null(parents)) {
      twin$intervene(target = parents, ifunction = condition[, parents, with = FALSE])
    }
    twin$simulate(n, fixedvars = twin$unames)
    for (j in 1:length(modellist)) {
      if( identical(modeltype[j],"predict")) {
        predictlist[[j]][, i] <- stats::predict(modellist[[j]], newdata = twin$simdata, ...)
      } else {
        predictlist[[j]][, i] <- do.call(modellist[[j]], list(newdata = twin$simdata, ...))
      }
    }
  }
  return(predictlist)
}



