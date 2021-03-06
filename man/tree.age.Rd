% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree.age.R
\name{tree.age}
\alias{tree.age}
\title{Calculating the age of nodes and tips in a tree.}
\usage{
tree.age(tree, age, order = "past")
}
\arguments{
\item{tree}{A \code{phylo} object.}

\item{age}{The age of the tree. If missing the age is set to be the tree height.}

\item{order}{Either "past" if the units express time since the present (e.g. million years ago), or "present" if the unit is expressed in time since the root.}
}
\description{
Calculates the age of each node and tip in a tree give the height of the tree or some specified age.
}
\examples{
## A dated random phylogeny with a root 50 units of time old.
tree.age(rtree(10), age = 50)
## A random tree with the distance since the root.
tree.age(rtree(10), order = 'present')

}
\seealso{
\code{\link{slice.tree}}, \code{\link{time.subsamples}}.
}
\author{
Thomas Guillerme
}
