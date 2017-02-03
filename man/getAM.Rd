\name{getAM}
\alias{getAM}
\title{
Construct Adjacency Matrix
}
\description{
Constructs a logical adjacency matrix from a SpatialPolygons* object. 
Polygons are only identified as adjacent if their intersection is a line or a polygon.
Hence, polygons that only touch at a single point are not considered adjacent.
}
\usage{
getAM(x)
}
\arguments{
  \item{x}{
A \code{SpatialPolygons*} object
}
}
\details{
}
\value{
A logical adjancency matrix.
}
\references{
}
\author{
Philipp Hunziker
}
\note{
}

\seealso{
}

\examples{
# Make chess board
library(sp)
gt <- GridTopology(c(0,0), c(1,1), c(8,8))
sg <- SpatialGrid(gt)
board <- as(as(sg, "SpatialPixels"), "SpatialPolygons")

# Get adjacency matrix
getAM(board)
}

\keyword{adjacency matrix}
