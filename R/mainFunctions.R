#' Title
#'
#' @param x
#' @param w
#' @param rho
#' @param lambda
#' @param H
#'
#' @return
#' @export
#'
#' @examples
Wfunc <- function(x,w,rho, lambda, H)
{
  if (x != w)
  {
    p1 <- (rho - w^2)/(rho - 2*w*x + x^2)
    p2 <- (w-x)^2/(rho - 2*w*x + x^2)
    x1 <- x
    x2 <- (rho-w*x)/(w-x)

    output <- p1*H(x1) + p2*H(x2)
  }

  if (x == w)
    output <- H(w) + lambda*(rho - w^2)

  return(output)
}

#' Title
#'
#' @param x1
#' @param w
#' @param rho
#' @param sig
#' @param nu
#'
#' @return
#' @export
#'
#' @examples
getOptimalDistribution = function(x1,w,rho)
{
  p1 <- (rho - w^2)/((rho - w^2) + (x1-w)^2)
  x2 <- w + (rho - w^2)/(w-x1)
  p2 <- 1 - p1

  output <- list(x = c(x1,x2), p = c(p1,p2))
  return(output)
}


#' Title
#'
#' @param H
#' @param nu
#' @param mu
#' @param sigma
#' @param lambda
#' @param direction
#'
#' @return
#' @export
#' @importFrom GenSA GenSA
#'
#' @examples
#' ####
#' #### Finding a the optimal probability measure (p,x) to the problem
#' ####
#' #### max P(X > c)
#' #### s.t. sum(p)  = 1
#' ####      sum(px) = 1
#' ####      sum(px^2) = 2
#' ####
#' #### where c is the 90th percentile of a standard exponential distribution
#' #### Note that the solution to this problem is known (see Theorem 3.3 of Bertsimas and Popescu)
#'
#' # Function and parameters for the integral of the objective
#' H <- function(x) as.numeric(qexp(0.9) <= x)
#' mu <- 1
#' sigma <- 2
#' lambda <- 0
#'
#' output <- computeBoundVal(H,mu, sigma, lambda)
#'
#' # Check that the bound matches the analytical solution in Bertsimas
#' CMsquare <- (sigma- mu^2)/mu^2
#' delta <-  qexp(0.9)/mu-1
#'
#' data.frame(Algorithm = output$bound, Analytical = CMsquare/(CMsquare + delta^2))
#'
#' # Check that the output is feasible
#' with(output$P, data.frame(moments = c(sum(p),sum(p*x), sum(p*x^2)), truth = c(1,mu,sigma) ))
computeBoundVal <- function(H,mu, sigma, lambda, direction = "max")
{
  scale <- if (direction == "min") 1 else -1

  if (mu^2 > sigma) output <- list(bound = scale*Inf, P = NULL)
  if (mu^2 == sigma)output <- list(bound = H(mu), P = list(p = 1, x = mu))
  if (mu^2 < sigma){
    Z = GenSAmodified(fn = function(x1) scale*Wfunc(x1, mu, sigma, lambda, H),
               par         = 0,
               lower       = 0,
               upper       = mu)

    output <- list(bound = scale*Z$value, P = getOptimalDistribution(x1 = Z$par,mu,sigma))
  }

  return(output)
}


#' Title
#'
#' @param H
#' @param nu
#' @param mu
#' @param sigma
#' @param lambda
#' @param direction
#'
#' @return
#' @export
#'
#' @examples
#' ####
#' #### Finding a the optimal probability measure (p,x) to the problem
#' ####
#' #### max P(X > c)
#' #### s.t. sum(p)  = 1
#' ####      sum(px) = 1
#' ####      sum(px^2) = 2
#' ####
#' #### where c is the 90th percentile of a standard exponential distribution
#' #### Note that the solution to this problem is known (see Theorem 3.3 of Bertsimas and Popescu)
#'
#' # Function and parameters for the integral of the objective
#' H <- function(x) as.numeric(qexp(0.9) <= x)
#' mu <- c(1,1+1e-5)
#' sigma <- c(2,2+1e-5)
#' lambda <- 0
#'
#' output <- computeBoundInt(H,mu, sigma, lambda)
#'
#' # Check that the bound is slightly larger (because we increase the feasible set) than the analytical solution in Bertsimas
#' CMsquare <- (sigma- mu^2)/mu^2
#' delta <-  qexp(0.9)/mu-1
#'
#' data.frame(Algorithm = output$bound, Analytical = CMsquare/(CMsquare + delta^2))
#'
#' # Check that the output is feasible
#' with(output$P, data.frame(lowerMomentBound = c(1, mu[1], sigma[1]), moments = c(sum(p),sum(p*x), sum(p*x^2)), upperMomentBound = c(1, mu[2], sigma[2]) ))
computeBoundInt <- function(H,mu, sigma, lambda, direction = "max")
{

  scale <- if (direction == "min") 1 else -1

  # Exclude trival scenarios
  if ( (mu[1] > mu[2]) || (sigma[1] > sigma[2])) output <- list(bound = scale*Inf, P = NULL)
  if (mu[1]^2 > sigma[2])  output <- list(bound = scale*Inf, P = NULL)
  if (mu[1]^2 == sigma[2]) output <- list(bound = H(mu), P = list(p = 1, x = mu))

  # Non-trivial scenario
  if (mu[1]^2 < sigma[2])
  {
    # First subprogram
    Z1 <- GenSAmodified(fn = function(args) scale*Wfunc(x = args[1], mu[2], rho = args[2], lambda, H),
               par         = c(0,max(sigma[1],mu[2]^2)),
               lower       = c(0,max(sigma[1],mu[2]^2)),
               upper       = c(mu[2],sigma[2])
               )
    P1 <- getOptimalDistribution(x1 = Z1$par[1],w = mu[2],rho = Z1$par[2])

    # Second Subprogram
    f2 <- function(args){ if (args[1] > args[2]) Inf else scale*Wfunc(x = args[1], w = args[2], sigma[2], lambda, H)}
    Z2 <- GenSAmodified(fn = f2,
               par   = c(0, mu[1]),
               lower = c(0, mu[1]),
               upper = rep(min(mu[2], sqrt(sigma[2])),2)
               )
    P2 <- getOptimalDistribution(x1 = Z2$par[1],w = Z2$par[2], rho =  sigma[2])

    bound <- max(scale*Z1$value,scale*Z2$value)
    P <- if (bound == scale*Z1$value) P1 else P2

    output <- list(bound = bound, P = P)
  }

  return(output)
}

#' Title
#'
#' @param H
#' @param paramOptim
#' @param direction
#'
#' @return
#' @export
#'
#' @examples
#' ###
#' ### Solve problem with equalities
#' ###
#' # Function and parameters for the integral of the objective
#' H <- function(x) as.numeric(qexp(0.9) <= x)
#' mu <- 1
#' sigma <- 2
#' lambda <- 0
#'
#' output <- computeBound(H,mu, sigma, lambda)
#'
#' # Check that the bound matches the analytical solution in Bertsimas
#' CMsquare <- (sigma- mu^2)/mu^2
#' delta <-  qexp(0.9)/mu-1
#'
#' data.frame(Algorithm = output$bound, Analytical = CMsquare/(CMsquare + delta^2))
#'
#' # Check that the output is feasible
#' with(output$P, data.frame(moments = c(sum(p),sum(p*x), sum(p*x^2)), truth = c(1,mu,sigma) ))
#'
#' ###
#' ### Solve problem with equalities
#' ###
#' mu <- c(1,1)
#' sigma <- c(2,2)
#' lambda <- 0
#'
#' output <- computeBound(H,mu, sigma, lambda)
#'
#' # Check that the bound is larger (because we increase the feasible set) than the analytical solution in Bertsimas
#' data.frame(Algorithm = output$bound, Analytical = CMsquare/(CMsquare + delta^2))
#'
#' # Check that the output is feasible
#' with(output$P, data.frame(lowerMomentBound = c(1, mu[1], sigma[1]), moments = c(sum(p),sum(p*x), sum(p*x^2)), upperMomentBound = c(1, mu[2], sigma[2]) ))
computeBound = function(H, mu, sigma, lambda,scale = 1, direction = "max")
{
  Nmu <- length(mu)
  Nsigma <- length(sigma)

  # Check parameters
  if (Nmu > 2 | Nmu < 1) return("mu must be either a scalar or a vector of length 2.")
  if (Nsigma > 2 | Nsigma < 1) return("sigma must be either a  scalar or vector of length 2.")

  # Compute bound
  if (Nmu == 1 & Nsigma == 1) output <- computeBoundVal(H, mu, sigma, lambda, direction)
  if (Nmu == 2 & Nsigma == 2) output <- computeBoundInt(H, mu, sigma, lambda, direction)
  if (Nmu == 1 & Nsigma == 2) output <- computeBoundInt(H, rep(mu,2), sigma, lambda, direction)
  if (Nmu == 2 & Nsigma == 1) output <- computeBoundInt(H, mu, rep(sigma,2), lambda, direction)

  output$bound <- scale*output$bound
  return(output)
 }
