#' Extension of the GenSA function available in the CRAN package of the same name
#'
#' This function is an extension of the function GenSA defined in the package GenSA available on CRAN.
#' The only difference with the existing version is that \emph{GenSAmodified} allows the user to pass
#' lower bounds and upper bounds of the optimization variable that are equal. Though this is a trivial
#' scenario in which case the optimal solution should be lower = upper, the current version of GenSA crashes.
#'
#' @inheritParams GenSA::GenSA
#' @importFrom GenSA GenSA
#' @examples
#'
#' # Test GenSAmodified in the univariate case
#' fn1 <- function(x) x
#' bound <- RobustTail:::GenSAmodified(fn = fn1,par = 1,lower = 1,upper = 1)
#' bound
#' GenSA::GenSA(fn = fn1,par = 1,lower = 1,upper = 1) # This crashes
#'
#' # Test GenSAmodified in the bivariate case with one lower bound equal to the upper bound
#' fn2 <- function(x) x[1] + x[2]
#' bound <- RobustTail:::GenSAmodified(fn = fn2,par = c(1,2),lower = c(1,1),upper = c(1,3))
#' bound$par
#' bound$value
#'
#' GenSA::GenSA(fn = fn2,par = c(1,2),lower = c(1,1),upper = c(1,3))# This crashes
#'
#' # Test GenSAmodified in the bivariate case with both lower bound are equal to the upper bounds
#' bound <- RobustTail:::GenSAmodified(fn = fn2,par = c(1,1),lower = c(1,1),upper = c(1,1))
#' bound$par
#' bound$value
#'
#' GenSA::GenSA(fn = fn2,par = c(1,1),lower = c(1,1),upper = c(1,1)) # This crashes
#'
#' # Let's check that when the lower bounds are strictly smaller than the upper bounds,
#' # all goes well
#' outGenSAmodified <- RobustTail:::GenSAmodified(fn = fn2,par = c(1,2),lower = c(1,2),upper = c(2,3))
#' outGenSA <- GenSA::GenSA(fn = fn2,par = c(1,2),lower = c(1,2),upper = c(2,3))
#'
#' # Check that they have the same optimal solution
#' all(outGenSAmodified$par == outGenSA$par)
#'
#' # Check that they have the same optimal objective value
#' outGenSAmodified$value == outGenSA$value

GenSAmodified <- function (par = NULL, fn, lower, upper, control = list(), ...)
{
  # Check necessary conditions to run GenSAmodified
  if (!is.function(fn) || is.null(fn)) {
    stop("'fn' has to be a R function.")
  }

  if (length(lower) != length(upper)) {
    stop("Lower and upper bounds vector do not have the same length")
  }

  # If all lower are different from the upper,
  # run the usual GenSA algorithm
  if (all(lower != upper)){
    output <- GenSA::GenSA(par = par, fn = fn, lower = lower, upper = upper, ...)
  }else {
    index <- which(lower == upper)
    newfn <- function(newx  = NULL) {
      x <- rep(NA,length(lower))
      x[index] <- lower[index]
      x[-index] <- newx
      output <- fn(x)
      return(output)
    }

    if (length(index) == length(lower)) {
      output <- list(value = newfn(), par = lower)
      return(output)
    }

    newlower <- lower[-index]
    newupper <- upper[-index]
    newpar <- par[-index]

    output <- GenSA::GenSA(par = newpar, fn = newfn, lower = newlower, upper = newupper, ...)

    # The output vector par of GenSA will be of the same dimension as the newlower bound
    # we need to include the value of lowerbound that was discarded in "par"
    if (is.list(output)){
      newpar <- rep(NA,length(lower))
      newpar[index] <- lower[index]
      newpar[-index] <- output$par

      output$par <- newpar
    }

  }

  return(output)
}
