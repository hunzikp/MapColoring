\name{getOptimalContrast}
\alias{getOptimalContrast}

\title{
Optimal Contrast Map Coloring
}
\description{
Chooses K colors from the \code{colors} vector and assigns them to the features in the map described by \code{x} such that (i) no two adjacent features have the same color and (ii) the average color contrast between adjacent features is maximized.
K is the minimal number of colors required to color the provided map in a manner such that no two adjacent features have the same color (i.e. an estimate of its chromatic number).
}
\usage{
getOptimalContrast(x, colors, beta = 1.5, seed = 1, ...)
}

\arguments{
  \item{x}{
An adjacency matrix or a \code{SpatialPolygons*} object
}
  \item{colors}{
A vector of color names. The function returns an error if the length of this vector is shorter than the number of colors required to color the map (K). Use \code{\link{getNColors}} to find K.
}
  \item{beta}{
A strictly positive parameter determining the degree to which low contrasts between adjacent colors are penalized. In most instances, the default value should work fine. 
}
  \item{seed}{
A random seed. Is required to make the results of the Cross-Entropy optimizer replicable.
}
  \item{\dots}{
Further arguments passed to \code{CEoptim}.
}
}
\details{
The number of colors required (K) and the assignment of non-adjacent color indices are determined using the DSATUR algorithm, see \code{\link{dsatur}}.
The assignment of actual colors from the \code{colors} vector to the color indices is determined by minimizing the sum of inverse contrast ratios between all adjacent features. 
Minimization is performed using the \code{CEoptim} Cross-Entropy optimizer, and the according package must be installed to use the \code{getOptimalContrast} function.
}
\value{
A vector of color names.
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

# Get candidate colors
colors <- rainbow(20)

# Optimal contrast color assignment
\dontrun{oc.colors <- getOptimalContrast(board, colors)}
\dontrun{plot(board, col=oc.colors)}
}

\keyword{map coloring}
