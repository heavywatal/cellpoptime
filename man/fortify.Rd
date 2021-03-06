% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fortify.R
\name{fortify_cellpop}
\alias{fortify_cellpop}
\alias{ggtree_fortify}
\alias{layout_rect}
\alias{layout_rect.list}
\alias{layout_rect.data.frame}
\title{Convert tbl_tree to plottable data.frame}
\usage{
fortify_cellpop(model, data, ..., method = "fdr", q = 0.05)

ggtree_fortify(data, ...)

layout_rect(data, ladderize = TRUE, ...)

\method{layout_rect}{list}(data, ladderize = TRUE, ...)

\method{layout_rect}{data.frame}(data, ladderize = TRUE, ...)
}
\arguments{
\item{model}{output of \code{scale_branches()}}

\item{data}{tbl_tree}

\item{...}{passed to \code{layout_rect} or \code{ggtree::fortify}.}

\item{method}{passed to \code{p.adjust()}}

\item{q}{threshold}

\item{ladderize}{Use \code{ape::ladderize} or not.}
}
\description{
Convert tbl_tree to plottable data.frame
}
\details{
\code{fortify_cellpop} prepares plottable data.frame.

\code{ggtree_fortify} is a wrapper of \code{ggtree::fortify}.

\code{layout_rect} is a simplified version of \code{ggtree::fortify},
whose output lacks "branch" and "angle" columns.
}
