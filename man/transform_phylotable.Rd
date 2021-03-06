% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trees.R
\name{transform_phylotable}
\alias{transform_phylotable}
\title{Transforms a table into a nwk string}
\usage{
transform_phylotable(phylotable)
}
\arguments{
\item{phylotable}{Phylo object or tibble from phylo. fmt as: parent,node,\<label\>}
}
\value{
A string with a rooted nwk string
}
\description{
This function takes in a phylo object or a tibble from a phylo object and
it returns its string representation in newick format.
The nodenumber is overwriten by the labels, if present.
}
