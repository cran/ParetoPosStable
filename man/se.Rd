\name{se}
\alias{se}
\alias{se.default}
%\alias{print.se}
\title{Approximated standard errors of Pareto Positive Stable (PPS) parameter estimates}
\description{
It approximates the stantard errors of PPS parameter estimates by bootstrapping.}
\usage{
  \method{se}{default}(PPSfit, k = 2000, show.iters = TRUE, ...)
}
\arguments{
  \item{PPSfit}{a \code{PPSfit} Object, typically from \code{PPS.fit()}.}
  \item{k}{the number of steps in the bootstrapping procedure.}
  \item{show.iters}{A logical argument specifying if the steps in the bootstrap iteration procedure are shown.}
  \item{\dots}{other arguments.}
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
se(fit, k = 50, show.iters = FALSE)
}
