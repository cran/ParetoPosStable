\name{dPPS}
\alias{dPPS}
\title{Density of the Pareto Positive Stable (PPS) distribution}
\description{\code{dPPS} returns values of the density for a PPS distribution}
\usage{
dPPS(x, lam, sc, v)
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{lam}{first shape parameter.}
  \item{sc}{scale parameter.}
  \item{v}{second shape parameter.}
}
\details{
See Sarabia and Prieto (2009) for the expression of the density.
}
\value{
The densities.
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
\code{\link{pPPS}} for the cumulative distribution function, \code{\link{qPPS}} for the quantile function, \code{\link{hPPS}} for the hazard function and \code{\link{rPPS}} for the random generation of data.
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
