\name{logLik.PPSfit}
\alias{logLik.PPSfit}
\title{Log-likelihood value of a PPSfit Object}
\description{It returns the log-likelihood value of a \code{PPSfit} Object}
\usage{
\method{logLik}{PPSfit}(object, \dots)
}
\arguments{
  \item{object}{a \code{PPSfit} Object.}
  \item{\dots}{other arguments.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{The log-likelihood.}
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
\code{\link{PPS.fit}}
}
\examples{
x <- rPPS(50, 1.2, 100, 2.3)
fit <- PPS.fit(x)
logLik(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
