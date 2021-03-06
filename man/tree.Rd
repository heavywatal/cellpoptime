% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree.R
\name{as_tbl_tree}
\alias{as_tbl_tree}
\alias{as_phylo}
\alias{as_multiphylo}
\alias{group_clade}
\alias{upstream_corner}
\title{Functions to process tidytree objects}
\usage{
as_tbl_tree(x)

as_phylo(x)

as_multiphylo(tbls)

group_clade(x, nodes)

upstream_corner(x, node)
}
\arguments{
\item{x}{data.frame}

\item{tbls}{list of tbl_tree}

\item{nodes}{integer IDs of ancestral nodes}

\item{node}{integer ID}
}
\description{
Functions to process tidytree objects
}
\details{
\code{as_tbl_tree} is a workaround until dplyr 0.8.0 is released.

\code{as_phylo} is a shortcut to convert data.frame to ape::phylo.

\code{as_multiphylo} adds class name.

\code{group_clade} is a simplified version of \code{tidytree::groupClade}.

\code{upstream_corner} finds the joint between a node and its parent
}
