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
#' @name R6causalimport
#' @references
#' J. Pearl (2009). Causality, 2nd edition, Cambridge University Press.
#' @import R6
#' @import data.table
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @importFrom stats dnorm
#' @importFrom stats predict
#' @importFrom MASS mvrnorm
#' @import igraph
#' @import dosearch
#' @import causaleffect
#' @import cfid
## @import qgraph
## @import sqldf
NULL


#' @include R6causal_utils.R
NULL

#' R6 Class for structural causal models
#'
#' An R6 class for structural causal models (SCM)
#' with latent variables and missing data mechanism. There are methods for
#' defining, printing, plotting, intervening and simulating SCMs.
#' @export
SCM <- R6::R6Class("SCM",
    lock_class = FALSE,
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
      .vstarnames = NULL,
      .rnames_prefix = NULL,
      .rnames_target = NULL,
      .uvnames = NULL,
      .vrnames = NULL,
      .uvrnames = NULL,
      .unames_dedicated = NULL,
      .unames_confounder = NULL,
      .dedicated_u = NULL,
      .is_linear_gaussian = FALSE,
      .rmapping_from_prefix = NULL,
      .vstarmapping_from_prefix = NULL,
      .rmapping_to_prefix = NULL,
      .vfsymb = NULL,
      .vrfsymb = NULL,
      .simdata = NULL,
      .simdata_obs = NULL,
      .adjmatrix = NULL,
      .igraph = NULL,
      .igraph_nodedicated = NULL,
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
        
        nchildren <- table( unlist( private$.vfsymb))
        usingle <- intersect(names(private$.uflist), names(nchildren)[nchildren == 1])
        private$.dedicated_u <- lapply(private$.vfsymb, FUN = intersect, y = usingle)
        private$.unames_dedicated <- unlist(private$.dedicated_u)
        private$.unames_confounder <- setdiff(private$.unames,private$.unames_dedicated)
        
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
        
        private$.igraph <- igraph::graph_from_adjacency_matrix( private$.adjmatrix, "directed")
        private$.igraph <- igraph::set.edge.attribute(private$.igraph, name = "description", value = "")
        private$.igraph <- igraph::set.edge.attribute(private$.igraph, name = "bidir_status", value = "")
        keepnodes <- igraph::V(private$.igraph)[union(private$.vnames, private$.unames_confounder)]
        private$.igraph_nodedicated <- igraph::induced.subgraph(private$.igraph, keepnodes)
        private$.igraph_bidirected <- private$.bidirected_presentation(private$.igraph_nodedicated,
                                                                        private$.unames_confounder) 
        private$.bidirected_edges_TF <- (igraph::edge.attributes( private$.igraph_bidirected)$description == "U")
        private$.bidirected_edges_TF_copy2 <-
          ( igraph::edge.attributes( private$.igraph_bidirected)$description == "U") &
          ( igraph::edge.attributes( private$.igraph_bidirected)$bidir_status == "copy2")
        private$.toporder <- names(igraph::topo_sort( private$.igraph, mode = "out"))
        private$.toporderv <- intersect( private$.toporder, private$.vnames)
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
        private$.toporder_md <- names(igraph::topo_sort( private$.igraph_md, mode = "out"))
        private$.topordervr_md <- intersect( private$.toporder_md, private$.vrnames)
        private$.toporderr_md <- intersect( private$.toporder_md, private$.rnames)

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
      .derive_SCM = function() {
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
          private$.rnames_target  <- names(private$.rflist)
          private$.rnames_prefix  <- paste0( private$.rprefix, names(private$.rflist))
          private$.rnames <- private$.rnames_prefix
          private$.rflist_prefix <- private$.rflist
          names(private$.rflist_prefix) <- private$.rnames
          private$.rmapping_from_prefix <- as.list(names(private$.rflist))
          names(private$.rmapping_from_prefix) <- private$.rnames
          private$.vstarmapping_from_prefix <- paste0(private$.rmapping_from_prefix, private$.starsuffix)
          names(private$.vstarmapping_from_prefix) <- private$.rnames
          private$.rmapping_to_prefix <- as.list(private$.rnames)
          names(private$.rmapping_to_prefix) <- names(private$.rflist)
          private$.vrflist <- c(private$.vflist, private$.rflist_prefix)
          private$.vrnames <- names(private$.vrflist)
          private$.uvrnames <- c(private$.unames, private$.vrnames)
          private$.derive_graph_md()
        }
      },
      .bidirected_presentation = function(G, l) {
        # modified from the package 'causaleffect'
        if(length(l) == 0) return(G)
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
            obs.new <- as.vector(t(as.matrix(expand.grid(inc.ind, out.ind))))
            n.obs.new <- length(obs.new)/2
            G <- G + igraph::edge(v[obs.new], description = rep(NA, n.obs.new)) # replace path v_1 -> L -> v_2 with v_1 -> v_2
          }
          if (out.len > 1) {
            unobs.new <- as.vector(combn(out.ind, 2))
            n.unobs.new <- length(unobs.new)/2
            #G <- G + edges(v[c(unobs.new, unobs.new[2:1, ])], description = rep("U", 2 * ncol(unobs.new))) # replace path v_1 <- L -> v_2 with v_1 <-> v_2
            G <- G + igraph::edge(v[unobs.new], description = rep("U", n.unobs.new),
                           bidir_status = rep("copy1", n.unobs.new))
            G <- G + igraph::edge(v[unobs.new[length(unobs.new):1]], description = rep("U", n.unobs.new),
                           bidir_status = rep("copy2", n.unobs.new))
          }
          if (unobs.len > 0 & out.len > 0) {
            unobs.old <- as.vector(t(as.matrix(expand.grid(unobs.ind, out.ind))))
            n.unobs.old <- length(unobs.old)/2
            #G <- G + edges(v[c(unobs.old, unobs.old[2:1, ])], description = rep("U", 2 * ncol(unobs.old))) # replace path v_1 <-> L -> v_2 with v_1 <-> v_2
            G <- G + igraph::edge(v[unobs.old], description = rep("U", n.unobs.old),
                           bidir_status = rep("copy1", n.unobs.old))
            G <- G + igraph::edge(v[unobs.old], description = rep("U", n.unobs.old),
                           bidir_status = rep("copy2", n.unobs.old))
          }
          G <- igraph::induced.subgraph(G, setdiff(v, l[i]))
        }
        return(G)
      },
      .parsefunction = function(obj) {
        if(is.function(obj)) {
          return(obj)
        }
        if(is.numeric(obj) | (is.data.frame(obj))) {
          return( constantfunction(as.numeric(obj)) )
        }
        if(is.factor(obj)) {
          return( constantfunction(obj) )
        }
        errormessage <- "A function must be specified as 1) a function or, 2) as numeric, 
        factor, or one-row data frame or, 3) as a character vector with the prespecified format 
        (see documentation)"
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
      #' @field vnames List of the names of observed variables.
      vnames = function(value) {
        if (missing(value)) {
          private$.vnames
        } else {
          stop("`$vnames` is read only", call. = FALSE)
        }
      },
      #' @field vstarnames List of the names of observed variables with NA's.
      vstarnames = function(value) {
        if (missing(value)) {
          private$.vstarnames
        } else {
          stop("`$vstarnames` is read only", call. = FALSE)
        }
      },
      #' @field vfsymb List of the arguments of structural functions of observed variables.
      vfsymb = function(value) {
        if (missing(value)) {
          private$.vfsymb
        } else {
          stop("`$vfsymb` is read only", call. = FALSE)
        }
      },
      #' @field uflist List of the structural functions of unobserved variables.
      uflist = function(value) {
        if (missing(value)) {
          private$.uflist
        } else {
          stop("`$uflist` is read only", call. = FALSE)
        }
      },
      #' @field unames List of the names of unobserved variables.
      unames = function(value) {
        if (missing(value)) {
          private$.unames
        } else {
          stop("`$unames` is read only", call. = FALSE)
        }
      },
      #' @field unames_dedicated List of the names of unobserved variables that have only one child.
      unames_dedicated = function(value) {
        if (missing(value)) {
          private$.unames_dedicated
        } else {
          stop("`$unames_dedicated` is read only", call. = FALSE)
        }
      },
      #' @field unames_confounder List of the names of unobserved variables that 
      #' have two or more children.
      unames_confounder = function(value) {
        if (missing(value)) {
          private$.unames_confounder
        } else {
          stop("`$unames_confounder` is read only", call. = FALSE)
        }
      },
      #' @field dedicated_u Named list of the names of unobserved variables that 
      #' have only one child which is the name of the element.
      dedicated_u = function(value) {
        if (missing(value)) {
          private$.dedicated_u
        } else {
          stop("`$dedicated_u` is read only", call. = FALSE)
        }
      },
      #' @field is_linear_gaussian Logical, does the SCM have linear functions 
      #' and Gaussian background variables?
      is_linear_gaussian = function(value) {
        if (missing(value)) {
          private$.is_linear_gaussian
        } else {
          # stopifnot(is.logical(value))
          # private$.is_linear_gaussian <- value
          # self
          stop("`$is_linear_gaussian` is read only", call. = FALSE)
        }
      },
      #' @field rflist List of the structural functions of missingness indicators.
      rflist = function(value) {
        if (missing(value)) {
          private$.rflist
        } else {
          stop("`$rflist` is read only", call. = FALSE)
        }
      },
      #' @field rfsymb List of the names of missingness indicators.
      rfsymb = function(value) {
        if (missing(value)) {
          private$.rfsymb
        } else {
          stop("`$rflist` is read only", call. = FALSE)
        }
      },
      #' @field rprefix Prefix used to mark missingness indicators.
      rprefix = function(value) {
        if (missing(value)) {
          private$.rprefix
        } else {
          stop("`$rprefix` is read only", call. = FALSE)
        }
      },
      #' @field starsuffix Suffix used to mark variables with missing data.
      starsuffix  = function(value) {
        if (missing(value)) {
          private$.starsuffix 
        } else {
          stop("`$starsuffix ` is read only", call. = FALSE)
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
      #' @field simdata_obs Data table containing data simulated from the SCM
      #' where missing values are indicated by \code{NA}.
      simdata_obs = function(value) {
        if (missing(value)) {
          private$.simdata_obs
        } else {
          stop("`$simdata_obs` is read only", call. = FALSE)
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
      #' @field igraph_nodedicated The graph of the SCM in the \code{igraph} form
      #' (without the dedicated U variables and the missing data mechanism).
      igraph_nodedicated = function(value) {
        if (missing(value)) {
          private$.igraph_nodedicated
        } else {
          stop("`$igraph_nodedicated` is read only", call. = FALSE)
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
     #' @param uflist A named list containing the functions for the background variables.
     #' @param vflist A named list containing the functions for the observed variables.
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
     # lg <- SCM$new("Linear Gaussian",linear_gaussian = list(nv = 50, edgeprob=0.5,
     #                vcoefdistr = function(n) {rnorm(n)},
     #                ccoefdistr = function(n) {rnorm(n)},
     #                ucoefdistr = function(n) {rnorm(n)}))
     initialize = function(name="An SCM", uflist = NULL, vflist = NULL, 
                           rflist = NULL, rprefix = "R_", starsuffix = "_md") {
        private$.name <- name
        private$.uflist <- lapply( uflist, private$.parsefunction)
        private$.vflist <- lapply( vflist, private$.parsefunction)
        if(!is.null(rflist)) {
          private$.rflist <- lapply( rflist, private$.parsefunction)
          private$.rprefix <- rprefix
          private$.starsuffix <- starsuffix
        }
        private$.derive_SCM()
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
#' @param subset Variable groups to be plotted: "uvr", "u2vr","vr","uv", "u2v" or "v".
#' @param ... other parameters passed to the plotting method
#' @examples
#' backdoor$plot()
#' backdoor$plot("v")
#SCM$set("public", "plot", 
plot = function(subset = "uvr", method = "igraph", ...) {
  if( !(method %in% c("qgraph","igraph"))) stop("The plotting method must be 'qgraph' or 'igraph'")
  if( !(subset %in% c("uvr","u2vr","vr","uv","u2v","v"))) {
    stop("The subset must be in c('uvr','u2vr','vr','uv','u2v','v')")
  }
  if( identical(subset,"uv") | identical(subset,"u2v") | 
      ((identical(subset,"uvr") | identical(subset,"u2vr")) & is.null(private$.rflist)))  {
    if(identical(subset,"uv") | identical(subset,"uvr")) {
      G <- private$.igraph
    } else {
      G <- igraph::induced.subgraph(private$.igraph, 
                                    union(private$.vnames,private$.unames_confounder))
    }
    if(identical(method,"igraph")) igraph::plot.igraph(G,...)
    if(identical(method,"qgraph")) {
      if (!requireNamespace("qgraph", quietly = TRUE)) {
        stop("Package \"qgraph\" needed for this plotting method to work. Please install it.",
             call. = FALSE)
      }
      qgraph::qgraph(
        input = igraph::as_edgelist(G, names = FALSE),
        labels = names(igraph::V(G)),
        edgelist = TRUE,
        weighted = FALSE,
        nNodes = length(igraph::V(G)),
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
  if( (identical(subset,"uvr") | identical(subset,"u2vr") | identical(subset,"vr")) & !is.null(private$.rflist) )  {
    if(identical(subset,"vr") | identical(subset,"uvr")) { 
      G <- private$.igraph_md
    } else {
      G <- igraph::induced.subgraph(private$.igraph_md, 
                                    union(private$.vnames,private$.unames_confounder))
    }
    if(identical(method,"igraph")) igraph::plot.igraph(G,...)
    if(identical(method,"qgraph")) {
      if (!requireNamespace("qgraph", quietly = TRUE)) {
        stop("Package \"qgraph\" needed for this plotting method to work. Please install it.",
             call. = FALSE)
      }
      qgraph::qgraph(
        input = igraph::as_edgelist(G, names = FALSE),
        labels = names(igraph::V(G)),
        edgelist = TRUE,
        weighted = FALSE,
        nNodes = length(igraph::V(G)),
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
    G <- igraph::induced.subgraph(private$.igraph_bidirected, private$.vnames)
    if(identical(method,"igraph")) igraph::plot.igraph(G,...)
    if(identical(method,"qgraph")) {
      if (!requireNamespace("qgraph", quietly = TRUE)) {
        stop("Package \"qgraph\" needed for this plotting method to work. Please install it.",
             call. = FALSE)
      }
      curvature <- 7.0 * private$.bidirected_edges_TF - 14.0 * private$.bidirected_edges_TF_copy2
      linetype <-  1 * (!private$.bidirected_edges_TF_copy2) * ( 1 + 1 * private$.bidirected_edges_TF)
      qgraph::qgraph(
        input = igraph::as_edgelist(G, names = FALSE),
        labels = names(igraph::V(G)),
        edgelist = TRUE,
        weighted = FALSE,
        nNodes = length(igraph::V(G)),
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
#Modified version of ChatGPT assisted translation of javascript written by Santtu Tikka
#' @description
#' Return a TikZ code for drawing the DAG of the SCM object in LaTeX.
#' @param subset Variable groups to be plotted: "uvr","vr","uv", or "v".
#' @param layoutfunction A layout function from igraph package.
#' @param labels A named list that gives the names of vertices in TikZ.
#' @param settings A list with the following elements:
#' @param ... Arguments to be passed to layoutfunction
#SCM$set("public", "tikz", 
tikz = function(subset = "uvr", 
                 layoutfunction = igraph::layout_with_lgl, 
                 labels = NULL,
                 settings = list(
                   force=FALSE,
                   borders=TRUE,
                   shape="circle",
                   size=5,
                   scale=2.0),
                 ...) {
          if( !(subset %in% c("uvr","vr","uv","u2v","v"))) stop("The subset must be in c('uvr','vr','uv','u2v','v')")
          if( (identical(subset,"uvr") | identical(subset,"vr")) & !is.null(private$.igraph_md) ) {
            G <- private$.igraph_md
          } else {
            if(identical(subset,"uv")) G <- private$.igraph
            if(identical(subset,"u2v")) G <-  igraph::induced.subgraph(private$.igraph, 
                                                                       union(private$.vnames,private$.unames_confounder))
            if(identical(subset,"v")) G <-  igraph::induced.subgraph(private$.igraph, 
                                                                     union(private$.vnames))
          }
          
          shape <- ifelse(settings$shape == "square", 
                          'regular polygon, regular polygon sides 4', 
                          settings$shape)
          style <- sprintf("{name = #1, %s, %s, inner sep = %fpt",
                           shape,
                           ifelse(settings$borders, "draw", ""),
                           settings$size)
          preamble <- "\\usepackage{tikz}\n\\usetikzlibrary{positioning, arrows.meta, shapes.geometric}\n" 
          preamble <- paste0(preamble, "\\tikzset{%\n  -Latex,semithick,\n  >={Latex[width=1.5mm,length=2mm]},\n")
          preamble <- ifelse(settings$force,
                             paste0(preamble, "  obs/.style n args = {2}", style, ", label = center:$#2$}\n"),
                             paste0(preamble, "  obs/.style = ", style, "}\n"))
          preamble <- paste0(preamble, "}")
          tikznodes <- ""
          tikzedges <- ""
          
          nodes <- names(igraph::V(G))
          if(is.null(labels)) {
            labels <- as.list(nodes)
            names(labels) <- nodes
          }
          edges <- igraph::as_edgelist(G)
          nodeCoordinates <- layoutfunction(G,...)
          if(is.list(nodeCoordinates)) nodeCoordinates <- nodeCoordinates$layout.dummy
          minY <- min(nodeCoordinates[,2])
          maxY <- max(nodeCoordinates[,2])
          minX <- min(nodeCoordinates[,1])
          maxX <- max(nodeCoordinates[,1])
          centerX <- 0.5 * (maxX + minX)
          centerY <- 0.5 * (maxY + minY)
          scaleX <- abs(maxX - minX)
          scaleY <- abs(maxY - minY)
          scale <- max(1, max(scaleX, scaleY))
          scale <- 12 / as.numeric(scale)
          tikznodes <- ""
          for (i in rev(seq_along(nodes))) {
            tikznodes <- paste(tikznodes, sprintf("\n  \\node [obs = {%s}%s] at (%f,%f) {$ %s $};",
                                                  nodes[i] ,
                                                  ifelse(settings$force, sprintf("{%s}", nodes[i] ), ""),
                                                  scale * (nodeCoordinates[i,1] - centerX),
                                                  -1 * scale * (nodeCoordinates[i,2] - centerY),
                                                  ifelse(settings$force, "\\vphantom{0}", labels[[ nodes[i] ]])))
          }
          tikzedges <- ""
          for (i in 1:nrow(edges)) {
            from <-  edges[i,1] 
            to <-  edges[i,2] 
            tikzedges <- paste(tikzedges, sprintf("\n  \\path [->] (%s) edge (%s);", from, to))
          }
          graph <-  paste("\\begin{tikzpicture}",tikznodes,tikzedges,"\n\\end{tikzpicture}")
          return(list(preamble = preamble, graph = graph))
},

#' @description
#' Return the parents of a set of vertices.
#' @param vnames A vector of vertex names
#' @param includeself Logical, should \code{vnames} to be included in the results (defaults TRUE)
#SCM$set("public", "pa", 
 pa = function(vnames, includeself = TRUE) {
          pa.ind <- unique(unlist(igraph::neighborhood(private$.igraph, 
                                                       order = 1, 
                                                       nodes = vnames, 
                                                       mode = "in")))
          pa <- V(private$.igraph)[pa.ind]$name
          if(!includeself) pa <- setdiff(pa, vnames)
          return(pa)
},
#' @description
#' Return the children of a set of vertices.
#' @param vnames A vector of vertex names
#' @param includeself Logical, should \code{vnames} to be included in the results (defaults TRUE)
#SCM$set("public", "ch", 
ch = function(vnames, includeself = TRUE) {
          ch.ind <- unique(unlist(igraph::neighborhood(private$.igraph, 
                                                       order = 1, 
                                                       nodes = vnames, 
                                                       mode = "out")))
          ch <- V(private$.igraph)[ch.ind]$name
          if(!includeself) ch <- setdiff(ch, vnames)
          return(de)
},
#' @description
#' Return the ancestors of a set of vertices.
#' @param vnames A vector of vertex names
#' @param includeself Logical, should \code{vnames} to be included in the results (defaults TRUE)
#SCM$set("public", "an", 
an = function(vnames, includeself = TRUE) {
          an.ind <- unique(unlist(igraph::neighborhood(private$.igraph, 
                                                       order = private$.num_uv, #vcount(private$.igraph), 
                                                       nodes = vnames, 
                                                       mode = "in")))
          an <- V(private$.igraph)[an.ind]$name
          if(!includeself) an <- setdiff(an, vnames)
          return(an)
},
#' @description
#' Return the descendants of a set of vertices.
#' @param vnames A vector of vertex names
#' @param includeself Logical, should \code{vnames} to be included in the results (defaults TRUE)
#SCM$set("public", "de", 
de = function(vnames, includeself = TRUE) {
          de.ind <- unique(unlist(igraph::neighborhood(private$.igraph, 
                                                       order = private$.num_uv, #vcount(private$.igraph), 
                                                       nodes = vnames, 
                                                       mode = "out")))
          de <- V(private$.igraph)[de.ind]$name
          if(!includeself) de <- setdiff(de, vnames)
          return(de)
},
#' @description
#' Add a new variable to the SCM object.
#' @param vfnew NULL or a named list containing the functions for the new observed variables.
#' @param ufnew NULL or a named list containing the functions for the new latent variables.
#' @param rfnew NULL or a named list containing the functions for the new missingness indicators.
#' @param rprefixnew NULL or the prefix of the missingness indicators.
#' @param starsuffixnew NULL orthe suffix for variables with missing data.
#' @examples
#' backdoor2 <- backdoor$clone()
#' backdoor2$add_variable(
#'    vfnew = list(
#'              w = function(uw, x) {
#'              return(as.numeric(uw < 0.4 + 0.3*x))}),
#'    ufnew = list(
#'             uw = function(n) {return(stats::runif(n))})
#' ) 
#SCM$set("public", "add_variable", 
add_variable = function(vfnew = NULL, ufnew = NULL, rfnew = NULL, 
                        rprefixnew = NULL,  starsuffixnew = NULL) {
  if( !is.null(vfnew)) {
    vf <- lapply( vfnew, private$.parsefunction)
    private$.vflist <- c(private$.vflist, vf)
  }
  if( !is.null(ufnew)) {
    uf <- lapply( ufnew, private$.parsefunction)
    private$.uflist <- c(private$.uflist, uf)
  }
  if( !is.null(rfnew)) {
    rf <- lapply( rfnew, private$.parsefunction)
    private$.rflist <- c(private$.rflist, rf)
  }
  if( !is.null(rprefixnew)) {
    private$.rprefix <- rprefixnew
  }
  if( !is.null(starsuffixnew)) {
    private$.starsuffix <- starsuffixnew
  }
  private$.derive_SCM()
},
#' @description
#' Remove variables from the SCM object.
#' @param variablenames Names of the variables to be removed.
#' @examples
#' backdoor2 <- backdoor$clone()
#' backdoor2$remove_variable(c("uy","y"))
#SCM$set("public", "remove_variable", 
remove_variable = function(variablenames) {
  if(length(setdiff(private$.vnames, variablenames))==0) { 
    stop("Cannot remove all observed variables.")
  }
  if( !is.null(private$.rflist)) {
    igraph <- private$.igraph_md
  } else {
    igraph <- private$.igraph
  }
  children_ind <- unlist(adjacent_vertices(igraph, 
                                           variablenames, mode = "out"))
  children <- setdiff(names(V(igraph))[children_ind], variablenames)
  if( length(children) > 0) { 
    stop("Cannot remove variables whose children are not removed.")
  }
  private$.vflist <- private$.vflist[setdiff(private$.vnames,variablenames)]
  private$.uflist <- private$.uflist[setdiff(private$.unames,variablenames)]
  if( !is.null(private$.rflist)) {
    private$.rflist <- private$.rflist[setdiff(private$.rnames,variablenames)]
    if(is.null(private$.rflist)) {
      private$.rprefix <- NULL
      private$.starsuffix <- NULL
      private$.rnames_target  <- NULL
      private$.rnames_prefix  <- NULL
      private$.rnames <- NULL
      private$.rflist_prefix <- NULL
      names(private$.rflist_prefix) <- NULL
      private$.rmapping_from_prefix <- NULL
      names(private$.rmapping_from_prefix) <- NULL
      private$.rmapping_to_prefix <- NULL
      names(private$.rmapping_to_prefix) <- NULL
      private$.vrflist <- NULL
      private$.vrnames <- NULL
      private$.uvrnames <- NULL
      private$.vrfsymb <- NULL
      private$.adjmatrix_md <- NULL
      private$.igraph_md <- NULL
      private$.toporder_md <- NULL
      private$.topordervr_md <- NULL
      private$.toporderr_md <- NULL
      private$.graphtext_md <- NULL
    }
  }
  private$.derive_SCM()
},

#' #' @include R6causal.R R6causal_examples.R
#' NULL
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
#SCM$set("public", "causal.effect", 
causal.effect = function(y,x,...) {
  return( causaleffect::causal.effect( y, x, G = private$.igraph_bidirected ,...))
},
#' @description
#' Is a causal effect or other query identifiable from given data sources?
#' Calls \code{\link[dosearch]{dosearch}} from the package \pkg{dosearch}.
#' See the documentation of \pkg{dosearch} for the details.
#' @param data Character string specifying the data sources.
#' @param query  Character string specifying the query of interest.
#' @param transportability Other parameters passed to \code{dosearch()}.
#' @param selection_bias Other parameters passed to \code{dosearch()}.
#' @param missing_data Other parameters passed to \code{dosearch()}.
#' @param control List of control parameters passed to \code{dosearch()}.
#' @return An object of class \code{dosearch::dosearch}.
#' @examples
#' backdoor$dosearch(data = "p(x,y,z)", query = "p(y|do(x))")
#SCM$set("public", "dosearch", 
dosearch = function(data, query, transportability = NULL, selection_bias = NULL,
                    missing_data  = NULL, control  = list()) {
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
},
#' @description
#' Is a counterfactual query identifiable from given data sources?
#' Calls \code{identifiable} from the package \pkg{cfid}.
#' See the documentation of \pkg{cfid} for the details.
#' @param gamma An R object that can be coerced into a \code{cfid::counterfactual_conjunction} object that represents the counterfactual causal query.
#' @param ... Other arguments passed to \code{cfid::identifiable}.
#' @return An object of class \code{cfid::query}.
#' @examples
#' backdoor$cfid(gamma = cfid::conj(cfid::cf("Y",0), cfid::cf("X",0, c(Z=1))) ) 
#SCM$set("public", "cfid", 
cfid = function(gamma,...) {
  if (!requireNamespace("cfid", quietly = TRUE)) {
    stop("Package \"cfid\" is needed to call method \"cfid\". Please install it.",
         call. = FALSE)
  }
  cfid_dag <- cfid::import_graph(private$.graphtext)
  return( cfid::identifiable(g = cfid_dag, gamma = gamma, ...))
},
     #' @description
     #' Apply an intervention to the SCM object.
     #' @param target Name(s) of the variables (in vflist, uflist or rflist) to be intervened.
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
        if(length(target) != length(ifunction)) {
            stop("The lengths of 'target' and 'ifunction' must be equal.")
        }
       if(length(ifunction) == 1) ifunction <- list(ifunction)
          for(i in 1:length(target)) {
            if(target[[i]] %in% private$.vnames) {
                private$.vflist[[ target[[i]]]] <- private$.parsefunction( ifunction[[i]])
            } else if(target[[i]] %in% private$.unames) {
                private$.uflist[[ target[[i]]]] <- private$.parsefunction( ifunction[[i]])
            } else if(target[[i]] %in% private$.rnames) {
                rtarget <- sub(paste0("^", private$.rprefix), "", target[[i]]) 
                private$.rflist[[ rtarget ]] <- private$.parsefunction( ifunction[[i]])
            } else {
              stop(paste("Unknown variable name:", target[[i]]))
            }
          }
        private$.derive_SCM()
      },
     #' @description
     #' Simulate data from the SCM object.
     #' Returns simulated data as a data.table and/or creates or updates \code{simdata} in the SCM object.
     #' If \code{no_missing_data = FALSE}, creates or updates also \code{simdata_obs}
     #' @param n Number of observations to be generated.
     #' @param no_missing_data Logical, should the generation of missing data skipped? (defaults FALSE).
     #' @param seed NULL or a number for \code{set.seed}.
     #' @param fixedvars List of variable names that remain unchanged or a data table/frame that contains the values of the fixed variables.
     #' @param store_simdata Logical, should the simulated data to be stored in the SCM object (defaults TRUE)
     #' @param return_simdata Logical, should the simulated data to be returned as the output (defaults FALSE)
     #' @examples
     #' backdoor$simulate(8, return_simdata = TRUE, store_simdata = FALSE)
     #' backdoor$simulate(10)
     #' backdoor$simdata
     simulate = function(n = 1, no_missing_data = FALSE, seed = NULL, fixedvars = NULL, store_simdata = TRUE, return_simdata = FALSE) {
        fixedvars_as_data <- is.data.frame(fixedvars) | is.data.table(fixedvars) 
        if(fixedvars_as_data) {
          fixedvarnames <- names(fixedvars)
        } else {
          fixedvarnames <- fixedvars
        }
        if(!is.null(seed)) set.seed(seed)
        if( is.null(fixedvars) | fixedvars_as_data) {
            private$.simdata <- data.table::data.table(matrix(as.numeric(NA), ncol = private$.num_uv, nrow = n))
            data.table::setnames(private$.simdata, private$.uvnames)
        } 
        for (i in 1:private$.num_u) {
            varchr <- private$.unames[i]
            if( !is.null(fixedvars)) {
              if(private$.unames[i] %in% fixedvarnames) { 
                if(!fixedvars_as_data) next
                if(fixedvars_as_data) {
                  data.table::set(private$.simdata, j = varchr, value = fixedvars[, ..varchr])
                } 
              } else {
                data.table::set(private$.simdata, j = varchr, value = private$.uflist[[i]](n = n))
              }
            } else {
              data.table::set(private$.simdata, j = varchr, value = private$.uflist[[i]](n = n))
            }
        }
        for (i in 1:private$.num_v) {
            varchr <- private$.toporderv[i]
            if( !is.null(fixedvars)) {
              if((varchr %in% fixedvarnames)) { 
                if(!fixedvars_as_data) next
                if(fixedvars_as_data) {
                  data.table::set(private$.simdata, j = varchr, value = fixedvars[, ..varchr])
                  next
                }
              }
            }
            arguments <- names(formals(private$.vflist[[ varchr ]]))
            if( identical(arguments,"...")) {
              testtry <- try(data.table::set(private$.simdata, j = varchr, value = private$.vflist[[ varchr ]]()))
            } else {
              testtry <- try(data.table::set(private$.simdata, j = varchr, value = do.call( private$.vflist[[ varchr ]],
                                                                            private$.simdata[ , ..arguments])))
            }
            if(inherits(testtry, "try-error")) {
              errormessage <- paste("Method simulate() failed when processing variable", 
                                  varchr,".", "The lower level error message was: ",testtry)
              stop(errormessage)
            }
          }
        data.table::setattr(private$.simdata, "SCMname", private$.name) 
        data.table::setattr(private$.simdata, "timestamp", date()) 
        data.table::setattr(private$.simdata, "seed", ifelse(is.null(seed),"",seed))
        if(no_missing_data | is.null(private$.rflist)) {
          data.table::setattr(private$.simdata, "graph", private$.graphtext)  
        } else {
          data.table::setattr(private$.simdata, "graph", private$.graphtext_md) 
        }
        if(return_simdata & (no_missing_data | is.null(private$.rflist))) {
          return(private$.simdata) 
        }
        if(!no_missing_data & !is.null(private$.rflist)) {
          private$.vstarnames <- paste0(private$.vnames, private$.starsuffix)
          private$.simdata[, c(private$.vstarnames)] <- private$.simdata[, private$.vnames, with = FALSE]
          for(i in 1:private$.num_r) {
            varchr <- private$.toporderr_md[i]
            if( !is.null(fixedvars)) {
              if((varchr %in% fixedvarnames)) {
                if(!fixedvars_as_data) next
                if(fixedvars_as_data) {
                  data.table::set(private$.simdata, j = varchr, value = fixedvars[, ..varchr])
                  next
                  }
                }
              }
              arguments <- names(formals(private$.rflist_prefix[[ varchr ]]))
              if( identical(arguments,"...")) {
                testtry <- try(data.table::set(private$.simdata, j = varchr, 
                                value = as.numeric(as.logical( private$.rflist_prefix[[ varchr ]]()))))
              } else {
                testtry <- try(data.table::set(private$.simdata, j = varchr, 
                                value = as.numeric(as.logical(do.call( private$.rflist_prefix[[ varchr ]],
                                private$.simdata[ , ..arguments, drop = FALSE])))))
              }
              if(inherits(testtry, "try-error")) {
                errormessage <- paste("Method simulate() failed when processing missingness indicator", 
                                      varchr,".", "The lower level error message was: ",testtry)
                stop(errormessage)
              }
            # TF_md <- as.logical( do.call( private$.rflist_prefix[[ varchr ]],
            #                               private$.simdata[ , ..arguments, drop=FALSE]))
            # data.table::set(private$.simdata, j = varchr, value = as.numeric(TF_md))
          }
          # Setting the value of an observed variable as NA if the missing indicator == 0
          #  which(! simdata_md[,j = ..varchr] ) tells the rows where this condition holds.
          for(i in 1:private$.num_r) {
            varchr <- private$.toporderr_md[i]
            if( !is.null(fixedvars)) {
              if((varchr %in% fixedvarnames)) next
            }
            data.table::set(private$.simdata, i = which(! private$.simdata[,j = ..varchr] ), 
                            j = private$.vstarmapping_from_prefix[[ varchr ]], value = NA)
          }
          private$.simdata_obs <- private$.simdata[, c(private$.vstarnames, private$.rnames), with = FALSE]
        }

        if(return_simdata & !no_missing_data & !is.null(private$.rflist)) {
          return(list( simdata = private$.simdata, simdata_obs = private$.simdata_obs)) 
        }
      }
    )
) #End of R6class SCM









