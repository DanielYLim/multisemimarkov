% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hazardLL.R
\name{hazardLL}
\alias{hazardLL}
\title{Log-Logistic hazard function}
\usage{
hazardLL(alpha1, beta1, t)
}
\arguments{
\item{alpha1}{a scale parameter.}

\item{beta1}{a shape parameter.}

\item{t}{time to event where t>0.}
}
\value{
hazard function of t with parameters alpha1 adn beta1.
}
\description{
Log-Logistic hazard function
}
\details{
Additional details...
 \deqn{h(t|alpha1, beta1) = ((beta1/alpha1)*(t/alpha1)^{beta1-1})/(1+(t/alpha1)^{beta1})}
}
\examples{
hazardLL(1, 2, 1.3)
\dontrun{
hazardLL(1, 2, 2.3)
}
}
