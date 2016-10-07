#' computeBoundVal
#'
#' This function returns the optimal objective value and solution of program (5) in the case
#' when nu := 1, and when the feasbile region is restricted to distribution functions with at most two point supports.
#'
#' @inheritParams W
#' @param mu :  the first order moment value in program (5)
#' @param sigma : the second order moment value in program (5)
#' @param direction : a string either "\emph{min}" or "\emph{max}" identifying the type of wether program (5)
#'  should be a min or a max program. Default is "max"
#'
#' @return  a list containing the optimal objective value \emph{bound} and distribution function \emph{P} of program (5). In particular, P is a list
#'  containing point p = (p1, p2) masses and supports points x =(x1,x2).
#' @export
#'
#' @examples
#' ####
#' #### We wish to solve
#' ####
#' #### max P(X > c)
#' #### s.t. sum(p)  = 1
#' ####      sum(px) = 1
#' ####      sum(px^2) = 2
#' ####
#' #### for some value c. We point out that the solution to this problem is known (see Theorem 3.3 of Bertsimas and Popescu),
#' ###
#'
#' # Function and parameters for the integral of the objective
#' c <- qexp(0.9)
#' H <- function(x) as.numeric(x >= c)
#' mu <- 1
#' sigma <- 2
#' lambda <- 0
#'
#' output <- computeBoundVal(H,mu, sigma, lambda)
#'
#' # Check that the bound matches the analytical solution in Bertsimas
#' CMsquare <- (sigma- mu^2)/mu^2
#' delta <- c/mu-1
#'
#' data.frame(Algorithm = output$bound, Analytical = CMsquare/(CMsquare + delta^2))
#'
#' # Check that the output is feasible
#' with(output$P, data.frame(moments = c(sum(p),sum(p*x), sum(p*x^2)), truth = c(1,mu,sigma) ))
computeBoundVal <- function(H,mu, sigma, lambda, direction = c("max", "min"))
{
  direction <- match.arg(direction)
  scale <- if (direction == "min") 1 else -1

  # Exclusion of trivial scenarios
  if (mu^2 > sigma) output <- list(bound = scale*Inf, P = NULL)
  if (mu^2 == sigma)output <- list(bound = H(mu), P = list(p = 1, x = mu))

  # Treatment of non-trivial scenarios
  if (mu^2 < sigma){
    Z <- GenSAmodified(fn = function(x1) scale*W(x1, mu, sigma, lambda, H),
                      par         = 0,
                      lower       = 0,
                      upper       = mu)
    output <- list(bound = scale*Z$value, P = getDistribution(x1 = Z$par,mu,sigma))
  }

  return(output)
}


#' computeBoundInt
#'
#' This function returns the optimal objective value and solution of program  (EC.19) in the case when nu := 1,
#' and when the feasbile region is restricted to distribution functions with at most two point supports.
#' @inheritParams computeBoundVal
#' @param H : the function defined in program  (EC.19)
#' @param mu : an ordered vector containing a lower and upper bound of the first order moment in program  (EC.19)
#' @param sigma : an ordered vector containing a lower and upper bound ofo the second order moment  in program  (EC.19)
#'
#' @return  a list containing the optimal objective value \emph{bound} and distribution function \emph{P} of program  (EC.19). In particular, P is a list
#'  containing point p = (p1, p2) masses and supports points x =(x1,x2).
#' @export
#'
#' @examples
#' ####
#' #### Finding a the optimal probability measure (p,x) to the problem
#' ####
#' #### max P(X > c)
#' #### s.t. sum(p)  = 1
#' ####      1 <= sum(px) <= 2
#' ####      2 <= sum(px^2) <= 3
#' ####
#' #### where c is some real number. The solution to the  problem
#' ####
#' #### max P(X > c)
#' #### s.t. sum(p)  = 1
#' ####      sum(px) = 1
#' ####      sum(px^2) = 2
#' #### is known (see Theorem 3.3 of Bertsimas and Popescu). Therefore, the optimal objective value
#' #### of the program with ineqality constraints should be an upper bound to the problem
#' #### with equality constraints.
#'
#' # Function and parameters for the integral of the objective
#' c <- qexp(0.9)
#' H <- function(x) as.numeric(c <= x)
#' mu <- c(1,2)
#' sigma <- c(2,3)
#' lambda <- 0
#'
#' output <- computeBoundInt(H,mu, sigma, lambda)
#'
#' # Check that the bound is larger (because we increase the feasible set) than the analytical solution in Bertsimas
#' sigma0 <- 2 ; mu0 <- 1
#' CMsquare <- (sigma0- mu0^2)/mu0^2
#' delta <-  c/mu0-1
#'
#' data.frame(Algorithm = output$bound, Analytical = CMsquare/(CMsquare + delta^2))
#'
#' # Check that the output is feasible
#' with(output$P, data.frame(lowerMomentBound = c(1, mu[1], sigma[1]), moments = c(sum(p),sum(p*x), sum(p*x^2)), upperMomentBound = c(1, mu[2], sigma[2]) ))
computeBoundInt <- function(H,mu, sigma, lambda, direction = c("max", "min"))
{
  direction <- match.arg(direction)
  scale <- if (direction == "min") 1 else -1

  # Exclude trival scenarios
  if ( (mu[1] > mu[2]) || (sigma[1] > sigma[2])) output <- list(bound = scale*Inf, P = NULL)
  if (mu[1]^2 > sigma[2])  output <- list(bound = scale*Inf, P = NULL)
  if (mu[1]^2 == sigma[2]) output <- list(bound = H(mu), P = list(p = 1, x = mu))

  # Non-trivial scenario
  if (mu[1]^2 < sigma[2])
  {
    # First subprogram
    lower <- c(0,max(sigma[1],mu[2]^2))
    upper <- c(mu[2],sigma[2])

    if (any(lower > upper)) { Z1 <- list(value = Inf) ; P1 <- list(NULL)
    } else{
      Z1 <- GenSAmodified(fn = function(args) scale*W(x = args[1], mu[2], rho = args[2], lambda, H),
                          par         = lower,
                          lower       = lower,
                          upper       = upper)

      P1 <- getDistribution(x1 = Z1$par[1],mu = mu[2],sigma = Z1$par[2])
    }


    # Second Subprogram
    lower <- c(0, mu[1])
    upper <- rep(min(mu[2], sqrt(sigma[2])),2)

    if (any(lower > upper)){ Z2 <- list(value = Inf) ; P2 <- list(NULL)
    } else{
      f2 <- function(args){ if (args[1] > args[2]) Inf else scale*W(x = args[1], w = args[2], sigma[2], lambda, H)}
      Z2 <- GenSAmodified(fn = f2,
                          par   = lower,
                          lower = lower,
                          upper = upper)
      P2 <- getDistribution(x1 = Z2$par[1],mu = Z2$par[2], sigma =  sigma[2])
    }

    bound <- max(scale*Z1$value,scale*Z2$value)
    P <- if (bound == scale*Z1$value) P1 else P2

    output <- list(bound = bound, P = P)
  }

  return(output)
}

#' computeBound
#'
#' This function is wrapper to solve either program (5) and (EC.19), in the case
#' when their feasbile region are restricted to distribution functions with at most two point supports.
#'
#' @inheritParams computeBoundVal
#' @param mu, sigma : sorted vectors of length at most 2 containing the lower and upper bound on the first moment in program (EC.19).
#' If both mu and sigma are scalars, then computeBound solves program (5), otherwise, computeBound solves program (EC.19).
#' @param scale : a scalar value by which we scale the optimal objective value. For instance, one can set scale := nu, where
#'  nu is as defined in program (5) and  (EC.19)
#' @param H : the function defined in the objective value of program  (5) and (EC.19)
#' @return a list containing the optimal objective value \emph{bound} and distribution function \emph{P} of program  (EC.19) or program (5). In particular, P is a list
#'  containing point p = (p1, p2) masses and supports points x =(x1,x2).
#' @export
#'
#' @examples
#' ####
#' #### We wish to solve
#' ####
#' #### max P(X > c)
#' #### s.t. sum(p)  = 1
#' ####      sum(px) = 1
#' ####      sum(px^2) = 2
#' ####
#' #### for some value c. We point out that the solution to this problem is known (see Theorem 3.3 of Bertsimas and Popescu).
#' #### We can use the function computeBound to solve this problem, which is of the form of program (5).
#'
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
#' ####
#' #### The previous problem is also equivalent to
#' ####
#' #### max P(X > c)
#' #### s.t. sum(p)  = 1
#' ####      1 <= sum(px) <= 1
#' ####      2 <= sum(px^2) <= 2
#' ####
#' #### We can use the function computeBound to solve this problem, which is of the form of program (EC.19).
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
#'
#' #### We now wish to solve
#' ####
#' #### Max P(x >b)
#' #### s.t. f(a) = eta
#' ####      f'(a) = -nu
#' ####      1-F(a) = beta
#' ####      f(x) is convex for all x => a
#' ####      f(x) is non-negative for all x => a
#' ####
#' #### This problem is of the form of program (1). By Theorem 4,  to solve that
#'
#'
#' a <- qexp(0.7)
#' eta <- dexp(a)
#' nu <- dexp(a)
#' beta <- 1-pexp(a)
#' mu <- eta/nu
#' sigma <- 2*beta/nu
#' lambda <- 1/2
#' b <- qexp(seq(0.7,0.99, by = 0.01))
#'
#' runFunc <- function(b){
#' H <- function(x) 1/2*(x - max(b-a,0) )^2*( max(b-a,0) <= x)
#' bound <- RobustTail::computeBound(H,mu,sigma,lambda,nu)$bound
#' output <- data.frame(b = b, bound = bound)
#' }
#'
#' dataPlot <- plyr::ldply(parallel::mclapply(X = b, FUN = runFunc, mc.cores = 4))
#'
#' library(ggplot2)
#' ggplot(dataPlot, aes(x = b, y = bound)) +
#' geom_line() +
#' ylim(c(0,0.3))
#'
#'
#'
#'
#'
#'
#'
computeBound = function(H, mu, sigma, lambda, nu = 1, direction = c("max", "min"))
{
  direction <- match.arg(direction)

  Nmu <- length(mu)
  Nsigma <- length(sigma)

  # Check parameters
  if (Nmu > 2 | Nmu < 1) return("mu must be either a scalar or a vector of length 2.")
  if (Nsigma > 2 | Nsigma < 1) return("sigma must be either a  scalar or vector of length 2.")

  # Check that nu >= 0
  scale <- if (direction == "min") 1 else -1
  if (nu < 0) return(list(bound = scale*Inf, P = NULL))

  # Compute bound
  if (Nmu == 1 & Nsigma == 1) output <- computeBoundVal(H, mu, sigma, lambda, direction)
  if (Nmu == 2 & Nsigma == 2) output <- computeBoundInt(H, mu, sigma, lambda, direction)
  if (Nmu == 1 & Nsigma == 2) output <- computeBoundInt(H, rep(mu,2), sigma, lambda, direction)
  if (Nmu == 2 & Nsigma == 1) output <- computeBoundInt(H, mu, rep(sigma,2), lambda, direction)

  output$bound <- nu*output$bound
  return(output)
}
