% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trees.R
\name{.printClade_wLength}
\alias{.printClade_wLength}
\title{Prints nwk format recursively}
\usage{
.printClade_wLength(phylotable, parent, depth = -1)
}
\arguments{
\item{phylotable}{Phylo object or tibble from phylo that is being queried}

\item{parent}{The node you want to investigate}

\item{depth}{OPTIONAL normally -1. This will return all. 0 returns nothing, one returns parent.}
}
\value{
A character string in nwk format
}
\description{
depth is the key! -1 to print everything, 0 prints nothing
}
