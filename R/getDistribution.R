#' getDistribution
#'
#' Returns a distribution function with two point masses
#' and supports (p1, p2, x1,x2) satisfying the three equations p1 + p2 = 1,  p1x1 + p2x2 = mu,
#' and m1x1^2 + p2x2^2 = sigma. This set of equations contains three unknown parameters (p1,p2, x2).
#' The parameters x1, mu, and sigma are given by the user.
#'
#' @param x1 : first point support
#' @param mu : mean of the distribution function
#' @param sigma : second moment of the distribution function
#'
#' @return a list containing the vectors p = (p1,p2) and x = (x1, x2)
#' @export
#'
#' @examples
#' x1 <- 1/2 ; mu <- 1 ; sigma <- 2
#' P <- getDistribution(x1,mu,sigma)
#' P
#'
#' # Check that the 3 equality constraints are satisfied
#' data.frame(moments = with(P,c(sum(p), sum(p*x), sum(p*x^2))), truth = c(1, mu,sigma))
#'
getDistribution = function(x1,mu,sigma)
{
  p1 <- (sigma - mu^2)/((sigma - mu^2) + (x1-mu)^2)
  x2 <- mu + (sigma - mu^2)/(mu-x1)
  p2 <- 1 - p1

  output <- list(x = c(x1,x2), p = c(p1,p2))
  return(output)
}

