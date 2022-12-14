% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{Cum_hazardUnif}
\alias{Cum_hazardUnif}
\title{Cumulative Uniform Hazard}
\usage{
Cum_hazardUnif(a14, b14, t)
}
\arguments{
\item{t}{}
}
\value{
cumulative hazard function of uniform distribution
}
\description{
Cumulative Uniform Hazard
}
\examples{
Cum_hazardUnif(1,4,2)
Cum_hazardUnif(1,4,0)  # outside of range gives 0
Cum_hazardUnif(1,4,5)  # outside of range gives 0
}
