\name{getColoring}
\alias{getColoring}

\title{
Return Map Coloring
}
\description{
Creates a map coloring for the provided object using the DSATUR graph coloring algorithm.
DSATUR assigns color indices in a manner such that no adjacent features have the same color, while trying to minimize the total number of colors used.
}
\usage{
getColoring(x)
}
\arguments{
  \item{x}{
An adjacency matrix or a \code{SpatialPolygons*} object.
}
}
\details{
In principle, the four-color-theorem ensures that the minimal number of colors required (or ``chromatic number'') is <5 for any planar graph. However, because DSATUR is greedy, the provided solution may not be globally optimal.
Note that the \code{getColoring} function only returns an integer vector representing color indices, not actual color values.}
\value{
An integer vector representing color indices.
}
\references{
}
\author{
Philipp Hunziker
}
\note{
}


\seealso{
\code{\link{getNColors}}
}

\examples{
# Make chess board
library(sp)
gt <- GridTopology(c(0,0), c(1,1), c(8,8))
sg <- SpatialGrid(gt)
board <- as(as(sg, "SpatialPixels"), "SpatialPolygons")

# Get color indices
getColoring(board)
}

