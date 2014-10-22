\name{pareto.fit}
\alias{pareto.fit}
\title{Fitting a Pareto distribution}
\description{It is an auxiliar function for fitting a Pareto distribution as a particular case of a Pareto Positive Stable distribution, allowing the scale parameter to be held fixed if desired.}
\usage{
pareto.fit(x, estim.method = "MLE", sigma = NULL, start, ...)
}
\arguments{
  \item{x}{the vector of observations.}
  \item{estim.method}{the estimation method, "MLE" or "OLS"}
  \item{sigma}{the value of the scale parameter, if it is known; if the value is \code{NULL}, the parameter is estimated.}
  \item{start}{unused argument from \code{PPS.fit}.}
  \item{\dots}{other arguments}
}
\details{This function is called by \code{PPS.fit()} when \code{Pareto} argument is \code{TRUE}.
}
\value{
A \code{PPSfit} Object.
}
\references{
Sarabia, J.M and Prieto, F. (2009). The Pareto-positive stable distribution: A new descriptive model for city size data, \emph{Physica A: Statistical Mechanics and its Applications}, \bold{388}(19), 4179-4191.
}

\seealso{
\code{\link{PPS.fit}}, \code{\link{coef.PPSfit}}, \code{\link{print.PPSfit}}, \code{\link{plot.PPSfit}}, \code{\link{GoF}}
}
