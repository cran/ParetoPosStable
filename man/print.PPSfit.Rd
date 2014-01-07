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
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
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
print(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
