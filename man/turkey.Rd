\name{turkey}
\alias{turkey}
\docType{data}
\title{
Population of Turkish cities and towns
}
\description{
Census population of Turkish cities and towns of more than 20,000 inhabitants in 1990 and 2000.
}
\usage{data(turkey)}
\format{
  A data frame with 280 observations on the following 4 variables.
  \describe{
    \item{\code{Name}}{cities and towns names.}
    \item{\code{Adm.}}{abbreviated province name of cities and towns.}
    \item{\code{Pop1990}}{the population in 1990.}
    \item{\code{Pop2000}}{the population in 2000.}
  }
}
\source{
\url{http://www.citypopulation.de/Turkey-C20.html}
}
\examples{
data(turkey)
summary(turkey)
}
\keyword{datasets}
