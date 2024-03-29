% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R6causal_linear.R
\name{analytic_linear_gaussian}
\alias{analytic_linear_gaussian}
\title{Simulate data from a conditional linear Gaussian SCM}
\usage{
analytic_linear_gaussian(scm, situation, n)
}
\arguments{
\item{scm}{A linear Gaussian SCM}

\item{situation}{A list with the following element:
\itemize{
\item condition : either a string that gives an SQL query ( e.g. "select x,y,z from DATA where" )
or a data.table consisting of the valid rows ( e.g. data.table::data.table( x = 0, y = 0))
}}

\item{n}{The number of rows in the data to be simulated}
}
\description{
Simulate data from a conditional linear Gaussian SCM
}
