\name{se.PPSfit}
\alias{se.PPSfit}
\title{Approximated standard errors of Pareto Positive Stable (PPS) parameter estimates}
\description{
It approximates the stantard errors of PPS parameter estimates by bootstrapping.}
\usage{
se.PPSfit(PPSfit, k = 2000, show.iters = TRUE)
}
\arguments{
  \item{PPSfit}{a \code{PPSfit} Object, typically from \code{PPS.fit()}.}
  \item{k}{the number of steps in the bootstrapping procedure.}
  \item{show.iters}{A logical argument specifying if the steps in the bootstrap iteration procedure are shown.}
}
\details{
The function simulates \code{k} samples from the model given in the \code{PPSfit} argument, fits them with the same method of estimation and uses the parameter estimates to approximate the standard errors.
}
\value{
A list with the standard errors.
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
coef(fit)
se.PPSfit(fit, k = 50, show.iters = FALSE)
}

