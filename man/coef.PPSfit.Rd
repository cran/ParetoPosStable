\name{coef.PPSfit}
\alias{coef.PPSfit}
\title{Parameter estimates in a PPSfit Object}
\description{\code{coef} returns the parameter estimates in a \code{PPSfit} Object}
\usage{
\method{coef}{PPSfit}(object, \dots)
}
\arguments{
  \item{object}{a \code{PPSfit} Object, typically from \code{PPS.fit()}.}
  \item{\dots}{other arguments.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{A vector vith the parameter estimates.}
%\references{
%% ~put references to the literature/web site here ~
%}
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
coef(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
