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

id <- function(x) x

cumsum0 <- function(x) {
  cs <- cumsum(x)
  return( c(0,cs[1:(length(cs)-1)]))
}

logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax <- function (x) {
  exp(x - logsumexp(x))
}

usample <- function(probs, categories, u) {
  lower <- t(apply(probs,1,cumsum0))
  upper <- t(apply(probs,1,cumsum))
  u_ind <- ( (u >= lower) & (u < upper))
  catmat <- matrix(categories, nrow(probs), ncol(probs), byrow = TRUE)
  as.vector(t(catmat))[as.vector(t(u_ind))]
}
#usample(matrix(c(0.2,0.3,0.5),10,3, byrow = TRUE), 1:3, runif(10))

rcateg <- function(x, u, softmax = FALSE, classnames = NULL, asfactor = TRUE, ordered = FALSE) {
  x <- rbind(x)
  if(softmax) {
    probs <- t(apply(x,1,softmax))
  } else {
    probs <- x
  }
  if(nrow(probs) == 1) probs <- cbind(rep(1,length(u))) %*% probs
  if(is.null(classnames)) classnames <- 1:ncol(x)
  simvar <- usample(probs, classnames, u)
  if(ordered) {
    return(as.ordered(simvar))
  }
  if(asfactor) {
    return(as.factor(simvar))
  } else {
    return(simvar) 
  }  
}
#rcateg(matrix(rnorm(30),10,3), runif(10))
#rcateg(matrix(rnorm(30),10,3), runif(10), ordered = TRUE)


# https://stackoverflow.com/questions/33850219/change-argument-names-inside-a-function-r
## Function to replace variables in function body
## expr is `body(f)`, keyvals is a lookup table for replacements
replace_vars <- function(expr, keyvals) {
  if (!length(expr)) return()
  for (i in seq_along(expr)) {
    if (is.call(expr[[i]])) expr[[i]][-1L] <- Recall(expr[[i]][-1L], keyvals)
    if (is.name(expr[[i]]) && deparse(expr[[i]]) %in% names(keyvals))
      expr[[i]] <- as.name(keyvals[[deparse(expr[[i]])]])
  }
  return( expr )
}

replace_formals <- function(forms, keyvals) {
  if (!length(forms)) return()
  formsnames <- names(forms)
  for (i in seq_along(forms)) {
    newname <- keyvals[formsnames[i]]
    if(!is.na(newname)) formsnames[i] <- newname
  }
  names(forms) <- formsnames
  return( forms )
}

rename_arguments <- function(f,newnames) {
  newbody <- replace_vars(body(f), newnames)
  newfunction <- f
  body(newfunction) <- newbody
  newformals <- replace_formals(formals(f), newnames) 
  formals(newfunction) <- newformals
  return(newfunction)
}
# f <- function(x, y) -x^2 + x + -y^2 + y
# newvals <- c('x'='x0', 'y'='y0')  # named lookup vector
# f2 <- rename_arguments(f,newvals)
# f3 <- rename_arguments(f,c('x'='x0'))