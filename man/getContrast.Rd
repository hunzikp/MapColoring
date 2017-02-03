\name{getContrast}
\alias{getContrast}

\title{
Color Contrast Ratio
}
\description{
Calculates the contrast ratio between colors given their relative luminance.
}
\usage{
getContrast(lum1, lum2)
}
\arguments{
  \item{lum1}{
Numeric vector of relative luminance values
}
  \item{lum2}{
Numeric vector of relative luminance values. Needs to be of same length as \code{lum1}
}
}
\details{
See http://www.w3.org/TR/WCAG20/#getContrast-ratiodef for definitions.
}
\value{
Numeric vector of contrast ratios.
}
\references{
\url{http://www.w3.org/TR/WCAG20/#getContrast-ratiodef}
}
\author{
Philipp Hunziker
}
\note{

}

\seealso{
\code{\link{getRelLum}}
}
\examples{
}
\keyword{contrast ratio}
