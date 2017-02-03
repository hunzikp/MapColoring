\name{MapColoring-package}
\alias{MapColoring-package}
\alias{MapColoring}
\docType{package}
\title{
\packageTitle{MapColoring}
}
\description{
\packageDescription{MapColoring}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{MapColoring}
\packageIndices{MapColoring}
}
\author{
\packageAuthor{MapColoring}

Maintainer: \packageMaintainer{MapColoring}
}
\references{
D.Brelaz (1979) - New Methods to color the vertices of a graph. Communications of the ACM: 22(4).
}

\keyword{ package }
\seealso{
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
