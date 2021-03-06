% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getDistribution.R
\name{getDistribution}
\alias{getDistribution}
\title{Returns a 2-point support feasible distribution function for program (5)}
\usage{
getDistribution(mu, sigma, x1 = NULL)
}
\arguments{
\item{mu}{First moment of the distribution function}

\item{sigma}{Second moment of the distribution function}

\item{x1}{A real number between 0 and mu. Default value is NULL, in which case,
x1 is drawn from a uniform distribution over [0,mu]}
}
\value{
a list composed of
\item{p}{a vector containing  the point masses }
\item{x}{a vector containing  the point support }
}
\description{
Returns a 2-point support feasible distribution function for program (5), i.e. a distribution function with
two point masses (p1,p2) and supports (x1,x2) satisfying the three equations p1 + p2 = 1,  p1x1 + p2x2 = mu,
and m1x1^2 + p2x2^2 = sigma.
}
\examples{
###
### Feasible system
###
mu <- 1 ; sigma <- 2
P <- getDistribution(mu,sigma, x1 = mu/2)
P

# Check that the 3 equality constraints are satisfied
data.frame(moments = with(P,c(sum(p), sum(p*x), sum(p*x^2))), truth = c(1, mu,sigma))

###
### Unfeasible system
###
 mu <- 2 ; sigma <- 2
getDistribution(mu,sigma)

###
### Unique solution
###
mu <- 1 ; sigma <- 1
getDistribution(mu,sigma)

}

