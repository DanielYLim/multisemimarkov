#
#' Simulation function for non-masked data
#'
#' @return
#' @export
#'
#' @examples
sim_no_masked <- function(n=400, alpha12 = 1, beta12 = 1, alpha13 = 2, beta13 = 2, alpha23 = 1, beta23 = 2, tau = 0.2, p = 0.7, C_max = 10, b=1){

  # data generation

  Data = data_generation_ver1(n=n, alpha12=alpha12, beta12=beta12, alpha13=alpha13, beta13=beta13, alpha23=alpha23, beta23=beta23,tau=tau, p=p, C_max=cmax)


  t1 = Data$t1
  t2 = Data$t2
  status12 = Data$status12
  status23 = Data$status23
  status13 = Data$status13



  cen12=1-sum(Data$status12)/n
  cen13=1-sum(Data$status13)/n
  cen23=1-sum(Data$status23)/n


  # estimation
  results<-nlminb(parms,log.lik.fun)
  hess = hessian(log.lik.fun, results$par)

  # variance estimates
  hess = solve(hess)
  se = sqrt(diag(hess))


  Out_data = list(data = Data, estimates = results, se_estimates = se)

  return(Out_data)
}
