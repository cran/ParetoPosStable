\name{qPPS}
\alias{qPPS}
\title{Quantile function of the the Pareto Positive Stable (PPS) distribution}
\description{\code{qPPS} returns the quantile function values for PPS distribution.}
\usage{
qPPS(p, lam, sc, v)
}
\arguments{
  \item{p}{vector of probabilities.}
  \item{lam}{first shape parameter.}
  \item{sc}{scale parameter.}
  \item{v}{second shape parameter.}
}
\details{
See Sarabia and Prieto (2009) for the expression of the quantile function.
}
\value{
The quantiles.
}
\references{
Sarabia, J.M and Prieto, F. (2009). The Pareto-positive stable distribution: A new descriptive model for city size data, \emph{Physica A: Statistical Mechanics and its Applications}, \bold{388}(19), 4179-4191.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dPPS}} for the density function, \code{\link{pPPS}} for the cumulative distribution function, \code{\link{hPPS}} for the hazard function, and \code{\link{rPPS}} for the random generation of data.
}
\examples{
print(x <- sort(rPPS(10, 1.2, 100, 2.3)))
dPPS(x, 1.2, 100, 2.3)
pPPS(x, 1.2, 100, 2.3)
qPPS(pPPS(x, 1.2, 100, 2.3), 1.2, 100, 2.3)
hPPS(x, 1.2, 100, 2.3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
