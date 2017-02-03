\name{dsatur}
\alias{dsatur}

\title{
Greedy DSATUR Graph Coloring Algorithm
}
\description{
R implementation of the greedy DSATUR graph coloring algorithm.
}
\usage{
dsatur(x, coloring=NULL)
}
\arguments{
  \item{x}{
Symmetric logical adjacency matrix.
}
  \item{coloring}{
An optional vector of integers of length \code{nrow(x)} indicating pre-colored vertices.
If provided, \code{dsatur} respects the initial coloring and just colors the remaining vertices.
Uncolored vertices are marked with a value of zero, whereas pre-colored vertices are identified with a positive integer.
Note that the set of pre-assigned colors should consist of an uninterrupted sequence of positive integers.
}
}
\details{
Partially ported from the Python implementation by Andrei Novikov (https://github.com/annoviko/pyclustering).
}
\value{
Integer vector indicating vertex colors.
}
\references{
D.Brelaz (1979) - New Methods to color the vertices of a graph. Communications of the ACM: 22(4).
}
\author{
Philipp Hunziker
}
\note{
}


\seealso{
}

\examples{
# Adjacency matrix for 2x2 grid
adj_mat <- matrix(c(FALSE, TRUE, TRUE, FALSE,
                    TRUE, FALSE, FALSE, TRUE,
                    TRUE, FALSE, FALSE, TRUE,
                    FALSE, TRUE, TRUE, FALSE), 4, 4)
# Get vertex coloring
dsatur(adj_mat)
}

\keyword{vertex coloring}
