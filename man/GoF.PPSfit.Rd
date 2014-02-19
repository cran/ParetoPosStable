\name{GoF.PPSfit}
\alias{GoF.PPSfit}
\title{Goodness of fit tests for the Pareto Positive Stable (PPS) distribution}
\description{Kolmogorov-Smirnov, Anderson-Darling and PPS goodness of fit tests to validate a PPS fit (typically from \code{PPS.fit()}).}
\usage{
GoF.PPSfit(PPSfit, k = 2000, show.iters = TRUE)
}
\arguments{
  \item{PPSfit}{a \code{PPSfit} Object.}
  \item{k}{the number of iterations in the bootstrap procedure to approximate the p-values.}
  \item{show.iters}{A logical argument specifying if the steps in the bootstrap iteration procedure are shown.}
}
\details{
It returns the Kolmogorov-Smirnov, the Anderson-Darling tests and a specific test for PPS distributions. p-values are approximated by a bootstrap procedure.

The specific goodness of fit test for PPS distributions is based on the linearity of the survival function vs. the scaled observations in a double log-log scale (see Sarabia and Prieto, 2009).
}
\value{
A list with the values of the tests statistics and the approximated p-values.
}
\references{
Sarabia, J.M and Prieto, F. (2009). The Pareto-positive stable distribution: A new descriptive model for city size data, \emph{Physica A: Statistical Mechanics and its Applications}, \bold{388}(19), 4179-4191.
}

\seealso{
\code{\link{PPS.fit}}, \code{\link{plot.PPSfit}}
}
\examples{
x <- rPPS(50, 1.2, 100, 2.3)
fit <- PPS.fit(x)
GoF.PPSfit(fit, k = 50, show.iters = FALSE)
}
