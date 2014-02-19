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

\value{A list with the parameter estimates.}

\seealso{
\code{\link{PPS.fit}}
}
\examples{
x <- rPPS(50, 1.2, 100, 2.3)
fit <- PPS.fit(x)
coef(fit)
}
