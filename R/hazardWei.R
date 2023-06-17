#' Weibull hazard function
#'
#' @param beta23 a shape parameter.
#' @param alpha23 a scale parameter.
#' @param t time to event where t>0.
#'
#' @return hazard function of t.
#' @export
#'
#' @examples
#' hazardWei(1,2,1.3)
#' \dontrun{
#' hazardWei(1,2,4)
#' }
hazardWei = function(beta23, alpha23, t){
  (beta23/alpha23)*(t/alpha23)**(beta23-1)
}
