#' Log-Logistic hazard function
#'
#' @param alpha1 a scale parameter.
#' @param beta1  a shape parameter.
#' @param t time to event where t>0.
#'
#' @return hazard function of t with parameters alpha1 adn beta1.
#' @export
#'
#' @details
#' Additional details...
#'  \deqn{h(t|alpha1, beta1) = ((beta1/alpha1)*(t/alpha1)^{beta1-1})/(1+(t/alpha1)^{beta1})}
#'
#'
#' @examples
#' hazardLL(1, 2, 1.3)
#' \dontrun{
#' hazardLL(1, 2, 2.3)
#' }
hazardLL <- function(alpha1, beta1, t) {
  ((beta1 / alpha1) * (t / alpha1)^(beta1 - 1)) / (1 + (t / alpha1)^beta1)
}
