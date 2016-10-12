#' Returns a 2-point support feasible distribution function for program (5)
#'
#' Returns a 2-point support feasible distribution function for program (5), i.e. a distribution function with
#' two point masses (p1,p2) and supports (x1,x2) satisfying the three equations p1 + p2 = 1,  p1x1 + p2x2 = mu,
#' and m1x1^2 + p2x2^2 = sigma.
#'
#' @param mu First moment of the distribution function
#' @param sigma Second moment of the distribution function
#' @param x1 A real number between 0 and mu. Default value is NULL, in which case,
#' x1 is drawn from a uniform distribution over [0,mu]
#'
#' @return a list composed of
#' \item{p}{a vector containing  the point masses }
#' \item{x}{a vector containing  the point support }
#' @export
#'
#' @examples
#' ###
#' ### Feasible system
#' ###
#' mu <- 1 ; sigma <- 2
#' P <- getDistribution(mu,sigma, x1 = mu/2)
#' P
#'
#' # Check that the 3 equality constraints are satisfied
#' data.frame(moments = with(P,c(sum(p), sum(p*x), sum(p*x^2))), truth = c(1, mu,sigma))
#'
#' ###
#' ### Unfeasible system
#' ###
#'  mu <- 2 ; sigma <- 2
#' getDistribution(mu,sigma)
#'
#' ###
#' ### Unique solution
#' ###
#' mu <- 1 ; sigma <- 1
#' getDistribution(mu,sigma)
#'
getDistribution = function(mu,sigma, x1 = NULL)
{
  if (mu^2 > sigma)  output <- NULL
  if (mu^2 == sigma) output <- list(p = 1, x = mu)
  if (mu^2 < sigma){

    if (is.null(x1) ) x1 <- runif(1,0,mu)

    p1 <- (sigma - mu^2)/((sigma - mu^2) + (x1-mu)^2)
    x2 <- mu + (sigma - mu^2)/(mu-x1)
    p2 <- 1 - p1

    output <- list(x = c(x1,x2), p = c(p1,p2))
  }
  return(output)
}

