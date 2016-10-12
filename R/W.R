#' Function W as defined in Theorem 5, in the case when nubar := 1.
#' @param x,w,rho Real numbers
#' @param lambda A real number equal to the limit of the ratio H(x)/x^2 as x goes to infinity
#' @param H The function H(x) as defined in program (5) and (EC.19)
#' @export
#' @examples
#' ###
#' ### Consider the case when :
#' ###    * h(x) := I(x >= c)
#' ###    * a : is the threshold a defined in program (1). We arbitrarily set a to
#' ###          the 70th percentile of a standard exponential
#' ###    * c : is arbitrarily fixed to the 90th percentile of a standard exponential
#' ###
#' ### We compute W(x1) as defined in Equation (8), i.e. in the case when
#' ###    * w := mu
#' ###    * rho := sigma
#' ###
#' a <- qexp(0.7)
#' c <- qexp(0.9)
#' H <- function(x) (x - max(c - a,0))^2/2*(x >= max(c - a,0))
#'
#' # Assume the true distribution function is a standard exponential
#' eta <- dexp(a)
#' nu <- dexp(a)
#' beta <- 1-pexp(a)
#'
#' mu <- eta/nu
#' sigma <- 2*beta/nu
#'
#' x1 <- seq(0,mu, length = 1000)
#' Wx1 <- sapply(x1, FUN = W, w = mu, rho = sigma, H = H, lambda = 1/2)
#' plot(x1,Wx1,type = "l", ylab = "W(x1)")
W <- function(x,w,rho, lambda, H)
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

