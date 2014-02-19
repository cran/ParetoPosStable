\name{print.PPSfit}
\alias{print.PPSfit}
\title{Printing a PPSfit Object}
\description{
It prints its argument (typically from \code{PPS.fit()}), returning some of the most important aspects.
}
\usage{
\method{print}{PPSfit}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{x}{a \code{PPSfit} object.}
  \item{digits}{the number of digits to be printed}
  \item{\dots}{other arguments.}
}

\references{
Sarabia, J.M and Prieto, F. (2009). The Pareto-positive stable distribution: A new descriptive model for city size data, \emph{Physica A: Statistical Mechanics and its Applications}, \bold{388}(19), 4179-4191.
}

\seealso{
\code{\link{PPS.fit}}
}
\examples{
x <- rPPS(50, 1.2, 100, 2.3)
fit <- PPS.fit(x)
print(fit)
}
