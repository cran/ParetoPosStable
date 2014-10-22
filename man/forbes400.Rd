\name{forbes400}
\alias{forbes400}
\docType{data}
\title{
Forbes 400 list. The Richest People in America in 2012
}
\description{
The Forbes 400 or 400 Richest Americans is a list published by Forbes magazine of the wealthiest 400 Americans, ranked by net worth. The dataset presents results about 2012.
}
\usage{data("forbes400")}
\format{
  A data frame with 400 observations on the following 2 variables.
  \describe{
    \item{\code{Name}}{members of the list names.}
    \item{\code{NetWorth}}{net worth in 2012 in $ billion.}
  }
}
\source{
\url{http://www.forbes.com/forbes-400/list/}
}
\examples{
data(forbes400)
summary(forbes400)
}
\keyword{datasets}
