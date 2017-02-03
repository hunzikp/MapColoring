\name{getNColors}
\alias{getNColors}
\title{
Get Required Number of Colors
}
\description{
Calculates the number of colors required to color the provided map, established using the DSATUR algorithm.
DSATUR assigns colors in a manner such that no adjacent features have the same color, while trying to minimize the total number of colors used.
}
\usage{
getNColors(x)
}
\arguments{
  \item{x}{
An adjacency matrix or a \code{SpatialPolygons*} object.
}
}
\details{
In principle, the four-color-theorem ensures that the minimal number of colors required (or ``chromatic number'') is <5 for any planar graph. However, because DSATUR is greedy, the provided solution may not be globally optimal in some instances.
}
\value{
A positive integer.
}
\references{
}
\author{
Philipp Hunziker
}
\note{
}


\seealso{
\code{\link{getColoring}}
}
\examples{
# Make chess board
library(sp)
gt <- GridTopology(c(0,0), c(1,1), c(8,8))
sg <- SpatialGrid(gt)
board <- as(as(sg, "SpatialPixels"), "SpatialPolygons")

# Get chromatic number
getNColors(board)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{chromatic number}