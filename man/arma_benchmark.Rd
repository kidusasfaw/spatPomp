% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arma_benchmark.R
\name{arma_benchmark}
\alias{arma_benchmark}
\title{Calculated log-ARMA log-likelihood benchmark for spatPomp models}
\usage{
arma_benchmark(spo, order = c(2, 0, 1))
}
\arguments{
\item{spo}{A spatPomp object}

\item{order}{A triple (p,d,q) for the ARIMA model fitted to the data. It is
intended that d=0}
}
\description{
Fits independent log-ARMA models for each unit, and calculates the conditional
log-likelihood for each observation, as well as log-likelihood for
each unit and total log-likelihood. A simple tool, but one with
practical applicability, as demonstrated by King et al (2008) and
Wheeler et al (2023). This function is designed for non-negative 
data, and adds 1 to each observation to avoid log(0).
}
\examples{
# Complete examples are provided in the package tests
\dontrun{
m <- he10(U = 5)
arma_benchmark(m)
}
}
\references{
\king2008

\wheeler2023
}
\seealso{
Other utilities: 
\code{\link{expand_params}()}
}
\author{
Edward L. Ionides
}
\concept{utilities}
