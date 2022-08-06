#' R6causal: R6 class for structural causal models
#'
#' Package R6causal implements an R6 class for structural causal models (SCM)
#' with latent variables and missing data mechanism.  The class contains methods for
#' 1) defining a structural causal model via functions, text or conditional probability tables,
#' 2) printing basic information on the model,
#' 3) plotting the graph for the model using packages `igraph` or `qgraph`,
#' 4) simulating data from the model, 5) applying an intervention,
#' 6) checking the identifiability of a query using the R packages `causaleffect` and `dosearch`,
#' 7) defining the missing data mechanism,
#' 8) simulating incomplete data from the model according to the specified missing data mechanism and
#' 9) checking the identifiability in a missing data problem using the R package `dosearch`.
#' In addition, there are functions for running experiments and doing counterfactual inference using simulation.
#'
#' @docType package
#' @name R6causal
#' @references
#' J. Pearl (2009). Causality, 2nd edition, Cambridge University Press.
#' @import R6
#' @import data.table
#' @importFrom stats runif
#' @import igraph
#' @import dosearch
#' @import causaleffect
## @import qgraph
## @import sqldf
NULL



# library(R6)
# library(igraph)
# library(qgraph)
# library(dosearch)
# library(causaleffect)
# library(sqldf)
# library(data.table)
# usethis::use_package("R6")
# usethis::use_package("data.table")
# usethis::use_package("stats")
# usethis::use_package("igraph")
# usethis::use_package("dosearch")
# usethis::use_package("causaleffect")
# #usethis::use_package("sqldf")
# usethis::use_package("sqldf","Suggests")
# usethis::use_package("qgraph","Suggests")
# #usethis::use_package("qgraph")
# usethis::use_vignette("using_R6causal")

# defaultW <- getOption("warn")
# options(warn = 4)


constantfunction <- function(constant) {
  force(constant)
  return( function(...) {  return(constant) } )
}

descendants <- function(node, G, includeself = TRUE) {
  # node  = a character vector of the names of the nodes
  # G = an igraph object
  de.ind <- unique(unlist(igraph::neighborhood(G, order = vcount(G), nodes = node, mode = "out")))
  de <- V(G)[de.ind]$name
  if(!includeself) de <- setdiff(de, node)
  return(de)
}

#' R6 Class for structural causal models
#'
#' An R6 class for structural causal models (SCM)
#' with latent variables and missing data mechanism. There are methods for
#' defining, printing, plotting, intervening and simulating SCMs.
#' @export
SCM <- R6::R6Class("SCM",
    private = list(
      .name = NULL,
      .uflist = NULL,
      .vflist = NULL,
      .rflist = NULL,
      .vrflist = NULL,
      .rprefix = NULL,
      .rflist_prefix = NULL,
      .starsuffix = NULL,
      .num_u = NULL,
      .num_v = NULL,
      .num_r = NULL,
      .num_uv = NULL,
      .num_vr = NULL,
      .num_uvr = NULL,
      .unames = NULL,
      .vnames = NULL,
      .rnames = NULL,
      .rnames_prefix = NULL,
      .rnames_target = NULL,
      .uvnames = NULL,
      .vrnames = NULL,
      .uvrnames = NULL,
      .rmapping_from_prefix = NULL,
      .rmapping_to_prefix = NULL,
      .vfsymb = NULL,
      .vrfsymb = NULL,
      .simdata = NULL,
      .simdata_md = NULL,
      .adjmatrix = NULL,
      .igraph = NULL,
      .igraph_bidirected = NULL,
      .bidirected_edges_TF = NULL,
      .bidirected_edges_TF_copy2 = NULL,
      .toporder = NULL,
      .toporderv = NULL,
      .graphtext = NULL,
      .adjmatrix_md = NULL,
      .igraph_md = NULL,
      .toporder_md = NULL,
      .toporderv_md = NULL,
      .toporderr_md = NULL,
      .topordervr_md = NULL,
      .graphtext_md = NULL,
      .graphtext_bidirected = NULL,
      .derive_graph = function() {
        private$.vfsymb <- vector( mode = "list", length(private$.vflist))
        names( private$.vfsymb) <- private$.vnames
        private$.adjmatrix <- matrix(0, private$.num_uv, private$.num_uv)
        colnames(private$.adjmatrix) <- private$.uvnames
        rownames(private$.adjmatrix) <- private$.uvnames
        for(i in 1:length(private$.vflist)) {
          varchr <- private$.vnames[i]
          arguments <- names(formals(private$.vflist[[i]]))
          private$.vfsymb[[i]] <- arguments
          if( !(identical(arguments,"..."))) {
            private$.adjmatrix[ arguments,  varchr] <- 1
          }
        }
        private$.igraph <- igraph::graph_from_adjacency_matrix( private$.adjmatrix, "directed")
        private$.igraph <- igraph::set.edge.attribute(private$.igraph, name = "description", value = "")
        private$.igraph <- igraph::set.edge.attribute(private$.igraph, name = "bidir_status", value = "")
        private$.igraph_bidirected <- private$.latent.projection(private$.igraph, private$.unames)
        private$.bidirected_edges_TF <- (igraph::edge.attributes( private$.igraph_bidirected)$description == "U")
        private$.bidirected_edges_TF_copy2 <-
          ( igraph::edge.attributes( private$.igraph_bidirected)$description == "U") &
          ( igraph::edge.attributes( private$.igraph_bidirected)$bidir_status == "copy2")
        private$.toporder <- igraph::topo_sort( private$.igraph, mode = "out")
        private$.toporderv <- intersect( names(private$.toporder), private$.vnames)

        nchildren <- table( unlist( private$.vfsymb))
        uarc <- intersect(names(private$.uflist), names(nchildren)[nchildren > 1])
        drawnvertex <- union( names( private$.vfsymb), uarc)
        graphtext <- ""
        for(i in 1: length( private$.vfsymb)) {
          for(j in 1:length( private$.vfsymb[[i]])) {
            if( private$.vfsymb[[i]][j] %in% drawnvertex ) {
              graphtext <- paste0( graphtext, private$.vfsymb[[i]][j], " -> ", private$.vnames[i], " \r\n ")
            }
          }
        }
        private$.graphtext <- graphtext
      },
      .derive_graph_md = function() {
        private$.vrfsymb <- vector( mode = "list", length(private$.vrflist))
        names( private$.vrfsymb) <- private$.vrnames
        private$.adjmatrix_md <- matrix(0, private$.num_uvr, private$.num_uvr)
        colnames(private$.adjmatrix_md) <- private$.uvrnames
        rownames(private$.adjmatrix_md) <- private$.uvrnames
        for(i in 1:length(private$.vrflist)) {
          varchr <- private$.vrnames[i]
          arguments <- names(formals(private$.vrflist[[i]]))
          private$.vrfsymb[[i]] <- arguments
          if( !(identical(arguments,"..."))) {
            private$.adjmatrix_md[ arguments,  varchr] <- 1
          }
        }
        private$.igraph_md <- igraph::graph_from_adjacency_matrix( private$.adjmatrix_md, "directed")
        private$.toporder_md <- igraph::topo_sort( private$.igraph_md, mode = "out")
        private$.topordervr_md <- intersect( names(private$.toporder_md), private$.vrnames)
        private$.toporderr_md <- intersect( names(private$.toporder_md), private$.rnames)

        nchildren <- table( unlist( private$.vrfsymb))
        uarc <- intersect(names(private$.uflist), names(nchildren)[nchildren > 1])
        drawnvertex <- union( names( private$.vrfsymb), uarc)
        graphtext <- ""
        for(i in 1: length( private$.vrfsymb)) {
          for(j in 1:length( private$.vrfsymb[[i]])) {
            if( private$.vrfsymb[[i]][j] %in% drawnvertex ) {
              graphtext <- paste0( graphtext, private$.vrfsymb[[i]][j], " -> ", private$.vrnames[i], " \r\n ")
            }
          }
        }
        private$.graphtext_md <- graphtext
      },
      .latent.projection = function(G, l) {
        # modified from the package 'causaleffect'
        to <- NULL
        from <- NULL
        description <- NULL
        for (i in 1:length(l)) {
          if( !("description" %in% names(igraph::edge_attr(G)))) stop("Edge attribute 'description' required for latent projection")
          e <- igraph::E(G)
          v <- igraph::get.vertex.attribute(G, "name")
          inc.edges <- e[.to(l[i]) & (is.na(description) | description != "U")]
          out.edges <- e[.from(l[i]) & (is.na(description) | description != "U")]
          unobs.edges <- e[.to(l[i]) & description == "U" & !is.na(description)]
          inc.ind <- igraph::get.edges(G, inc.edges)[ ,1]
          out.ind <- igraph::get.edges(G, out.edges)[ ,2]
          unobs.ind <- setdiff( igraph::get.edges(G, unobs.edges)[ ,1], out.ind)
          inc.len <- length(inc.ind)
          out.len <- length(out.ind)
          unobs.len <- length(unobs.ind)
          if (inc.len > 0 & out.len > 0) {
            obs.new <- t(as.matrix(expand.grid(inc.ind, out.ind)))
            G <- G + igraph::edge(v[c(obs.new)], description = rep(NA, ncol(obs.new))) # replace path v_1 -> L -> v_2 with v_1 -> v_2
          }
          if (out.len > 1) {
            unobs.new <- combn(out.ind, 2)
            #G <- G + edges(v[c(unobs.new, unobs.new[2:1, ])], description = rep("U", 2 * ncol(unobs.new))) # replace path v_1 <- L -> v_2 with v_1 <-> v_2
            G <- G + igraph::edge(v[c(unobs.new[1,], unobs.new[2, ])], description = rep("U", ncol(unobs.new)),
                           bidir_status = rep("copy1", ncol(unobs.new)))
            G <- G + igraph::edge(v[c(unobs.new[2,], unobs.new[1, ])], description = rep("U", ncol(unobs.new)),
                           bidir_status = rep("copy2", ncol(unobs.new)))
          }
          if (unobs.len > 0 & out.len > 0) {
            unobs.old <- t(as.matrix(expand.grid(unobs.ind, out.ind)))
            #G <- G + edges(v[c(unobs.old, unobs.old[2:1, ])], description = rep("U", 2 * ncol(unobs.old))) # replace path v_1 <-> L -> v_2 with v_1 <-> v_2
            G <- G + igraph::edge(v[c(unobs.old[1,], unobs.old[2, ])], description = rep("U", ncol(unobs.old)),
                           bidir_status = rep("copy1", ncol(unobs.old)))
            G <- G + igraph::edge(v[c(unobs.old[2,], unobs.old[1, ])], description = rep("U", ncol(unobs.old)),
                           bidir_status = rep("copy2", ncol(unobs.old)))
          }
          G <- igraph::induced.subgraph(G, setdiff(v, l[i]))
          e.dat <- as.data.frame(igraph::get.edges(G, igraph::E(G)))
          e.dat[ ,3:4] <- igraph::edge.attributes(G)
          G <- igraph::subgraph.edges(G, which(!duplicated(e.dat)), delete.vertices = FALSE)
        }
        return(G)
      },
      .parsefunction = function(obj) {
        if(is.function(obj)) {
          return(obj)
        }
        if(is.numeric(obj)) {
          return( constantfunction(obj) )
        }
        errormessage <- "A function must be specified as a function or as numeric or as a character vector with the prespecified format (see documentation)"
        if(is.character(obj)) {
          textlist <- unlist( strsplit(obj, ":",  fixed = TRUE))
          if( length(textlist) != 2) stop( errormessage )
          functiontext <- paste0("function(", textlist[1],") { return(", textlist[2],") }")
          return( eval( parse( text = functiontext))) # Default envir = parent.frame . Is this correct?
        } else {
          stop( errormessage )
        }
      }
    ),
    active = list(
      #' @field vflist List of the structural functions of observed variables.
      vflist = function(value) {
        if (missing(value)) {
          private$.vflist
        } else {
          stop("`$vflist` is read only", call. = FALSE)
        }
      },
      #' @field vfsymb List of the names of observed variables.
      vfsymb = function(value) {
        if (missing(value)) {
          private$.vfsymb
        } else {
          stop("`$vflist` is read only", call. = FALSE)
        }
      },
      #' @field simdata Data table containing data simulated from the SCM.
      simdata = function(value) {
        if (missing(value)) {
          private$.simdata
        } else {
          stopifnot( length(union(names(value), names(private$.simdata))) ==
                       length(intersect(names(value), names(private$.simdata)))   )
          private$.simdata <- value
          self
        }
      },
      #' @field simdata_md Data table containing data simulated from the SCM
      #' where missing values are indicated by \code{NA}.
      simdata_md = function(value) {
        if (missing(value)) {
          private$.simdata_md
        } else {
          stop("`$simdata_md` is read only", call. = FALSE)
        }
      },
      #' @field igraph The graph of the SCM in the \code{igraph} form
      #' (without the missing data mechanism).
      igraph = function(value) {
        if (missing(value)) {
          private$.igraph
        } else {
          stop("`$igraph` is read only", call. = FALSE)
        }
      },
      #' @field igraph_bidirected The graph of the SCM in the \code{igraph} form where
      #' latent variables are presented by bidirected arcs.
      igraph_bidirected = function(value) {
        if (missing(value)) {
          private$.igraph_bidirected
        } else {
          stop("`$igraph_birected` is read only", call. = FALSE)
        }
      },
      #' @field igraph_md The graph of the SCM in the \code{igraph} form including
      #' the missing data mechanism.
      igraph_md = function(value) {
        if (missing(value)) {
          private$.igraph_md
        } else {
          stop("`$igraph_md` is read only", call. = FALSE)
        }
      },
      #' @field toporder A vector giving the topological order of variables.
      toporder = function(value) {
        if (missing(value)) {
          private$.toporder
        } else {
          stop("`$toporder` is read only", call. = FALSE) #TODO: check that the given order is topological
        }
      },
      #' @field toporderv A vector giving the topological order of observed
      #' variables.
      toporderv = function(value) {
        if (missing(value)) {
          private$.toporderv
        } else {
          stop("`$toporderv` is read only", call. = FALSE)
        }
      },
      #' @field graphtext A character string that gives the edges of the graph
      #' of the SCM (without the missing data mechanism).
      graphtext = function(value) {
        if (missing(value)) {
          private$.graphtext
        } else {
          stop("`$graphtext` is read only", call. = FALSE)
        }
      },
      #' @field graphtext_md A character string that gives the edges of the graph
      #' of the SCM including the missing data mechanism.
      graphtext_md = function(value) {
        if (missing(value)) {
          private$.graphtext_md
        } else {
          stop("`$graphtext_md` is read only", call. = FALSE)
        }
      },
      #' @field name The name of the SCM.
      name = function(value) {
        if (missing(value)) {
          private$.name
        } else {
          stopifnot(is.character(value), length(value) == 1)
          private$.name <- value
          self
        }
      }
    ),
   public = list(
     #' @description
     #' Create a new SCM object.
     #' @param name Name.
     #' @param uflist A named list containing the functions for latent variables.
     #' @param vflist A named list containing the functions for observed variables.
     #' @param rflist A named list containing the functions for missingness indicators.
     #' @param rprefix The prefix of the missingness indicators.
     #' @param starsuffix The suffix for variables with missing data.
     #' @return A new `SCM` object.
     #' @examples
     #' backdoor <- SCM$new("backdoor",
     #'  uflist = list(
     #'   uz = function(n) {return(stats::runif(n))},
     #'   ux = function(n) {return(stats::runif(n))},
     #'   uy = function(n) {return(stats::runif(n))}
     #'  ),
     #'  vflist = list(
     #'   z = function(uz) {
     #'     return(as.numeric(uz < 0.4))},
     #'   x = function(ux, z) {
     #'     return(as.numeric(ux < 0.2 + 0.5*z))},
     #'   y = function(uy, z, x) {
     #'     return(as.numeric(uy < 0.1 + 0.4*z + 0.4*x))}
     #'  )
     #' )
     initialize = function(name, uflist, vflist, rflist = NULL, rprefix = "R_", starsuffix = "_md") {
        private$.name <- name
        private$.uflist <- lapply( uflist, private$.parsefunction)
        private$.vflist <- lapply( vflist, private$.parsefunction)
        if(!is.null(rflist)) {
          private$.rflist <- lapply( rflist, private$.parsefunction)
        }
        private$.num_u <- length(uflist)
        private$.num_v <- length(vflist)
        private$.num_r <- length(rflist)
        private$.num_uv <- private$.num_u + private$.num_v
        private$.num_vr <- private$.num_v + private$.num_r
        private$.num_uvr <- private$.num_u + private$.num_v + private$.num_r
        private$.vnames <- names(private$.vflist)
        private$.unames <- names(private$.uflist)
        private$.uvnames <- c( private$.unames, private$.vnames)
        private$.derive_graph()
        if( !is.null(rflist)) {
          private$.rprefix <- rprefix
          private$.starsuffix <- starsuffix
          private$.rnames_target  <- names(private$.rflist)
          private$.rnames_prefix  <- paste0( rprefix, names(private$.rflist))
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
      },
     #' @description
     #' Print a summmary of the SCM object.
     #' @examples
     #' backdoor
     print = function() {
        cat("Name of the model: ", private$.name,"\n\n")
        cat("Graph: \n")
        cat("",private$.graphtext,"\n")
        cat("Functions of background (exogenous) variables: \n\n")
        print(private$.uflist)
        cat("Functions of endogenous variables: \n\n")
        print(private$.vflist)
        cat("Topological order of endogenous variables: \n")
        print(private$.toporderv)
        cat("\n")
        if(is.null(private$.rflist)) {
          cat("No missing data mechanism")
        } else {
          cat("Functions of missingness indicators: \n\n")
          print(private$.rflist)
        }
      },
     #' @description
     #' Plot the DAG of the SCM object.
     #' @param method Plotting method: "qgraph" or "igraph".
     #' @param subset Variable groups to be plotted: "uvr","vr","uv", or "v".
     #' @param ... other parameters passed to the plotting method
     #' @examples
     #' backdoor$plot()
     #' backdoor$plot("v")
      plot = function(subset = "uvr", method = "igraph", ...) {
        if( !(method %in% c("qgraph","igraph"))) stop("The plotting method must be 'qgraph' or 'igraph'")
        if( !(subset %in% c("uvr","vr","uv","v"))) stop("The subset must be in c('uvr','vr','uv','v')")
        if( identical(subset,"uv") | (identical(subset,"uvr") & is.null(private$.rflist)))  {
          if(identical(method,"igraph")) igraph::plot.igraph(private$.igraph,...)
          if(identical(method,"qgraph")) {
            if (!requireNamespace("qgraph", quietly = TRUE)) {
              stop("Package \"qgraph\" needed for this plotting method to work. Please install it.",
                   call. = FALSE)
            }
            qgraph::qgraph(
              input = igraph::as_edgelist(private$.igraph, names = FALSE),
              labels = names(igraph::V(private$.igraph)),
              edgelist = TRUE,
              weighted = FALSE,
              nNodes = length(igraph::V(private$.igraph)),
              rescale = TRUE,
              aspect = TRUE,
              normalize = TRUE,
              curveAll = TRUE,
              curveScale = TRUE,
              curveScaleNodeCorrection = FALSE,
              edge.color = "black",
              ...
            )
          }
        }
        if( (identical(subset,"uvr") | identical(subset,"vr")) & !is.null(private$.rflist) )  {
          if(identical(method,"igraph")) igraph::plot.igraph(private$.igraph_md,...)
          if(identical(method,"qgraph")) {
            if (!requireNamespace("qgraph", quietly = TRUE)) {
              stop("Package \"qgraph\" needed for this plotting method to work. Please install it.",
                   call. = FALSE)
            }
            qgraph::qgraph(
              input = igraph::as_edgelist(private$.igraph_md, names = FALSE),
              labels = names(igraph::V(private$.igraph_md)),
              edgelist = TRUE,
              weighted = FALSE,
              nNodes = length(igraph::V(private$.igraph_md)),
              rescale = TRUE,
              aspect = TRUE,
              normalize = TRUE,
              curveAll = TRUE,
              curveScale = TRUE,
              curveScaleNodeCorrection = FALSE,
              edge.color = "black",
              ...
            )
          }
        }
        if( identical(subset,"v"))  {
          if(identical(method,"igraph")) igraph::plot.igraph(private$.igraph_bidirected,...)
          if(identical(method,"qgraph")) {
            if (!requireNamespace("qgraph", quietly = TRUE)) {
              stop("Package \"qgraph\" needed for this plotting method to work. Please install it.",
                   call. = FALSE)
            }
            curvature <- 7.0 * private$.bidirected_edges_TF - 14.0 * private$.bidirected_edges_TF_copy2
            linetype <-  1 * (!private$.bidirected_edges_TF_copy2) * ( 1 + 1 * private$.bidirected_edges_TF)
          qgraph::qgraph(
            input = igraph::as_edgelist(private$.igraph_bidirected, names = FALSE),
            labels = names(igraph::V(private$.igraph_bidirected)),
            edgelist = TRUE,
            weighted = FALSE,
            nNodes = length(igraph::V(private$.igraph_bidirected)),
            rescale = TRUE,
            aspect = TRUE,
            normalize = FALSE,
            curveAll = TRUE,
            curveScale = TRUE,
            curveScaleNodeCorrection = FALSE,
            curve = curvature,
            lty = linetype,
            edge.color = "black",
            ...
          )
          }
        }
      },
     #' @description
     #' Apply an intervention to the SCM object.
     #' @param target Name(s) of the variables in vflist to be intervened.
     #' @param ifunction Either numeric value(s) or new structural function(s) for the target variables.
     #' @examples
     #' # A simple intervention
     #' backdoor_x1 <- backdoor$clone()  # making a copy
     #' backdoor_x1$intervene("x",1) # applying the intervention
     #' backdoor_x1$plot() # to see that arrows incoming to x are cut
     #'
     #' # An intervention that redefines a structural equation
     #' backdoor_yz <- backdoor$clone()  # making a copy
     #' backdoor_yz$intervene("y",
     #'     function(uy, z) {return(as.numeric(uy < 0.1 + 0.8*z ))}) # making y a function of z only
     #' backdoor_yz$plot() # to see that arrow x -> y is cut
     intervene = function(target, ifunction) {
        #Only interventions to v variables are currently implemented
        if(length(target) != length(ifunction)) {
            stop("The lengths of 'target' and 'ifunction' must be equal.")
        }
        if( length(target) > 1) {
          for(i in 1:length(target)) {
            private$.vflist[[ target[[i]]]] <- private$.parsefunction( ifunction[[i]])
          }
        } else {
          private$.vflist[[ target ]] <- private$.parsefunction(ifunction)
        }
        private$.derive_graph()
      },
     #' @description
     #' Simulate data from the SCM object.
     #' Creates or updates \code{simdata}.
     #' If \code{no_missing_data = FALSE}, creates or updates also \code{simdata_md}
     #' @param n Number of observations to be generated.
     #' @param no_missing_data Logical, should the generation of missing data skipped? (defaults FALSE).
     #' @param fixedvars List of variables that remain unchanged.
     #' @examples
     #' backdoor$simulate(10)
     #' backdoor$simdata
      simulate = function(n = 1, no_missing_data = FALSE, fixedvars = NULL) {
          if( is.null(fixedvars)) {
            simdata <- data.table::data.table(matrix(as.numeric(NA), ncol = private$.num_uv, nrow = n))
            colnames(simdata) <- private$.uvnames
          } else {
            simdata <- private$.simdata
          }
          for (i in 1:private$.num_u) {
            if( !is.null(fixedvars)) {
              if(private$.unames[i] %in% fixedvars) next
            }
            set(simdata, j = i, value = private$.uflist[[i]](n = n))
          }
          for (i in 1:private$.num_v) {
            varchr <- private$.toporderv[i]
            if( !is.null(fixedvars)) {
              if((varchr %in% fixedvars)) next
            }
            arguments <- names(formals(private$.vflist[[ varchr ]]))
            if( identical(arguments,"...")) {
              set(simdata, j = varchr, value = private$.vflist[[ varchr ]]())
            } else {
              set(simdata, j = varchr, value = do.call( private$.vflist[[ varchr ]],
                                                        simdata[ , ..arguments]))
            }
          }
          private$.simdata <- simdata
          private$.simdata_md <- NULL

        if(!no_missing_data & !is.null(private$.rflist)) {
          if( is.null(fixedvars) | !any(private$.rnames %in% fixedvars)) {
            simdata_md <- cbind( private$.simdata, data.table::data.table(matrix(NA, ncol = private$.num_r, nrow = n)))
            colnames(simdata_md)  <-  c( colnames(private$.simdata), private$.rnames)
          } else {
            simdata_md <- private$.simdata_md
          }

          for(i in 1:private$.num_r) {
            varchr <- private$.toporderr_md[i]
            if( !is.null(fixedvars)) {
              if((varchr %in% fixedvars)) next
            }
            arguments <- names(formals(private$.rflist_prefix[[ varchr ]]))
            TF_md <- as.logical( do.call( private$.rflist_prefix[[ varchr ]],
                                         simdata_md[ , ..arguments, drop=FALSE]))
            set(simdata_md, j = varchr, value = as.numeric(TF_md))
            #set(simdata_md, i = which(!TF_md), j = private$.rmapping_from_prefix[[ varchr ]], value = NA)
          }
          for(i in 1:private$.num_r) {
            varchr <- private$.toporderr_md[i]
            if( !is.null(fixedvars)) {
              if((varchr %in% fixedvars)) next
            }
            set(simdata_md, i = which(! simdata_md[,j = ..varchr] ) , j = private$.rmapping_from_prefix[[ varchr ]], value = NA)
          }

          mdvars <- c(private$.rnames_target, private$.rnames_prefix)
          simdata_md <- simdata_md[, ..mdvars, drop = FALSE]
          names(simdata_md) <- c(paste0(private$.rnames_target, private$.starsuffix), private$.rnames_prefix)
          private$.simdata_md <- simdata_md
          }
      },
     #' @description
     #' Is a causal effect identifiable from observational data?
     #' Calls the implementation of ID algorithm from package \pkg{causaleffect}.
     #' See the documentation of \code{\link[causaleffect]{causal.effect}} for the details.
     #' @param y A vector of character strings specifying target variable(s).
     #' @param x A vector of character strings specifying intervention variable(s).
     #' @param ... Other parameters passed to \code{\link[causaleffect]{causal.effect}}.
     #' @return An expression for the joint distribution of the set of variables (y) given
     #' the intervention on the set of variables (x) conditional on (z) if the effect is
     #'  identifiable. Otherwise an error is thrown describing the graphical structure
     #'  that witnesses non-identifiability.
     #'  @examples
     #'  backdoor$causal.effect(y = "y", x = "x")
      causal.effect = function(y,x,...) {
        return( causaleffect::causal.effect( y, x, G = private$.igraph_bidirected ,...))
      },
     #' @description
     #' Is a causal effect or other query identifiable from given data sources?
     #' Calls \code{\link[dosearch]{dosearch}} from the package \pkg{dosearch}.
     #' See the documentation of dosearch for the details.
     #' @param data Character string specifying the data sources.
     #' @param query  Character string specifying the query of interest.
     #' @param transportability Other parameters passed to \code{dosearch()}.
     #' @param selection_bias Other parameters passed to \code{dosearch()}.
     #' @param missing_data Other parameters passed to \code{dosearch()}.
     #' @param control List of control parameters passed to \code{dosearch()}.
     #' @return An object of class \code{dosearch}.
     #' @examples
     #' backdoor$dosearch(data = "p(x,y,z)", query = "p(y|do(x))")
      dosearch = function(data, query, transportability = NULL, selection_bias = NULL,
                          missing_data  = NULL, control  = NULL) {
        if(is.null(missing_data)  & !is.null(private$.rflist)) {
          missing_data <- paste( paste0(private$.rnames_prefix," : ",private$.rnames_target), collapse = ", ")
        }

        if(is.null(missing_data)) {
          graph <- private$.graphtext
        } else {
          graph <- private$.graphtext_md
        }
        return( dosearch::dosearch(data = data, query = query, graph = graph,
                                     transportability = transportability, selection_bias = selection_bias,
                                     missing_data = missing_data, control = control))
      }
    )
)


cumsum0 <- function(x) {
 cs <- cumsum(x)
 return( c(0,cs[1:(length(cs)-1)]))
}

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

#' SCM "backdoor" used in the examples.
#'
#' Variable z fulfills the back-door criterion for P(y|do(x))
#' @examples
#' backdoor
#' backdoor$plot()
#' @export
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

#' SCM "frontdoor" used in the examples.
#'
#' Variable z fulfills the front-door criterion for P(y|do(x))
#' @examples
#' frontdoor
#' frontdoor$plot()
#' @export
frontdoor <- SCM$new("frontdoor",
                     uflist = list(
                       uz = function(n) {return(stats::runif(n))},
                       ux = function(n) {return(stats::runif(n))},
                       uy = function(n) {return(stats::runif(n))},
                       u = function(n) {return(stats::runif(n))}
                     ),
                     vflist = list(
                       x = function(ux,u) {
                         return( as.numeric(ux < 0.2 + 0.7*u))},
                       z = function(uz,x) {
                         return( as.numeric(uz < 0.1 + 0.7*x))},
                       y = function(uy,z,u) {
                         return( as.numeric(uy < 0.1 + 0.3*z + 0.5*u + 0.2*u*z))}
                     )
)

#' SCM "trapdoor" used in the examples.
#'
#' Variable z is a trapdoor variable for P(y|do(x))
#' @references
#' J. Helske, S. Tikka, J. Karvanen (2021). Estimation of causal effects with
#' small data in the presence of trapdoor variables,
#' Journal of the Royal Statistical Society Series A, 184(3), 1030-1051,
#' http://doi.org/10.1111/rssa.12699
#' @examples
#' trapdoor
#' trapdoor$plot()
#' @export
trapdoor <- SCM$new("trapdoor",
                    uflist = list(
                      uw = function(n) {return(stats::runif(n))},
                      uz = function(n) {return(stats::runif(n))},
                      ux = function(n) {return(stats::runif(n))},
                      uy = function(n) {return(stats::runif(n))},
                      uwx = function(n) {return(stats::runif(n))},
                      uwy = function(n) {return(stats::runif(n))}
                    ),
                    vflist = list(
                      w = function(uw,uwx,uwy) {
                        return( as.numeric(uw < 0.3 + 0.3*uwx + 0.3*uwy))},
                      z = function(uz,w) {
                        return( as.numeric(uz < 0.4 + 0.4*w))},
                      x = function(ux,z,uwx) {
                        return( as.numeric(ux < 0.3 + 0.5*z - 0.2*uwx ))},
                      y = function(uy,x,uwy) {
                        return( as.numeric(uy < 0.1 + 0.4*x + 0.4*uwy))}
                    )
)

#' SCM "backdoor_md" used in the examples.
#'
#' Variable z fulfills the back-door criterion for P(y|do(x)).
#' Variable z is missing completely at random. The missingness of
#' variables x and y depend on z.
#' @examples
#' backdoor_md
#' backdoor_md$plot()
#' @export
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


#' Counterfactual inference via simulation
#'
#' @param scm An SCM object
#' @param situation A list or a character string. The list has the following elements:
#' \itemize{
#' \item do : NULL or a list containing named elements 'target' and 'ifunction' that
#' specify the intervention carried out in the situation
#' \item condition : either a string that gives an SQL query ( "select x,y,z from DATA where" )
#' or a data.table consisting of the valid rows
#' }
#' The character string specifies an SQL query ( "select x,y,z from DATA where" )
#' @param target A vector of variable names that specify the target
#' variable(s) of the counterfactual intervention.
#' @param ifunction A list of functions for the counterfactual intervention.
#' @param n Size of the data to be simulated
#' @param slotsize A scalar, the number of rows to be simulated for a slot
#' @param maxslots A scalar, the maximum number of slots
#' @param situationSQL Logical, is the situation defined as an SQL query, defaults FALSE
#' @return A data table representing the situation after the counterfactual intervention
#' @examples
#' cfdata <- counterfactual(backdoor,
#'                          situation = list(do = list(target = "x", ifunction = 0),
#'                          condition = data.table::data.table( x = 0, y = 0)),
#'                          target = "x",
#'                          ifunction = 1,
#'                          n = 100000)
#' mean(cfdata$y)
#' @export
counterfactual <- function(scm, situation, target, ifunction, n, slotsize = 10000, maxslots = 100, situationSQL = FALSE) {
  twin <- scm$clone()
  if(!is.null( situation$do)) {
    twin$intervene(situation$do$target, situation$do$ifunction)
  }
  twin$simulate(1)
  cf_data <- twin$simdata[0]
  sloti <- 1
  while( nrow(cf_data) < n & sloti <= maxslots) {
    twin$simulate(slotsize)
    if(situationSQL) {
      if (!requireNamespace("sqldf", quietly = TRUE)) {
        stop("Package \"sqldf\" needed when situationSQL == TRUE. Please install it.",
             call. = FALSE)
      }
      DATA <- twin$simdata
      validdata <- sqldf::sqldf(situation$condition)
    } else {
      validdata <- merge( twin$simdata, situation$condition, by = names(situation$condition), allow.cartesian = TRUE)
    }
    cf_data <- rbind(cf_data, validdata)
    sloti <- sloti + 1
  }
  truen <- min(n, nrow(cf_data))
  twin$simdata <- cf_data[1:truen, ]
  twin$intervene(target, ifunction)
  fixedvars <- names( twin$simdata)
  fixedvars <- setdiff( fixedvars, descendants(target, scm$igraph) )
  twin$simulate(truen, fixedvars =  fixedvars)
  return( twin$simdata)
}

#' Conduct a sequence of interventions and collect the simulated data.
#'
#' @param SCM An SCM object
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
run_experiment <- function(SCM, intervene, response, n) {
  response_list <- vector(mode = "list", length = length(response))
  names(response_list) <- response
  icombs <- do.call(CJ,intervene)
  resptemp <- data.table::as.data.table( matrix(NA, nrow = n, ncol = nrow(icombs)))
  for (j in 1:length(response_list)) {
     response_list[[j]] <- resptemp
  }
  for (i in 1:nrow(icombs)) {
    SCMcopy <- SCM$clone()
    SCMcopy$intervene( names(intervene), as.numeric(icombs[i,]))
    SCMcopy$simulate(n)
    for (j in 1:length(response)) {
      response_list[[ response[j] ]][,i] <- SCMcopy$simdata[, get(response[j])]
    }
  }
  return(list( interventions = icombs, response_list = response_list))
}


