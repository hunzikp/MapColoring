\name{getRelLum}
\alias{getRelLum}
\title{
Calculate Relative Luminocity
}
\description{
Calculates the relative luminocity for a matrix of colors specified in sRGB.
}
\usage{
getRelLum(color, maxColorValue = 255)
}
\arguments{
  \item{color}{
An sRGB matrix of colors
}
  \item{maxColorValue}{
Maximum value on the color bands (usually 255 or 1)
}
}
\details{
See http://www.w3.org/TR/WCAG20/#relativeluminancedef for definitions.
}
\value{
A numeric vector of relative luminocity values.
}
\references{
\url{http://www.w3.org/TR/WCAG20/#relativeluminancedef}
}
\author{
Philipp Hunziker
}
\note{
}

\seealso{
}
\examples{
}
\keyword{ luminocity }
