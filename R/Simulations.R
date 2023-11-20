# Simulation

#' Estimation unmasked data
#'
#' @param B
#' @param n
#' @param alpha12
#' @param beta12
#' @param alpha13
#' @param beta13
#' @param alpha23
#' @param beta23
#' @param tau
#' @param p
#' @param C_max
#' @param alpha04
#' @param beta04
#' @return nlm output with some information of generated data
#' @export
#' @examples
#' B=10
#' n=200
#' alpha12 = 2
#' beta12 = 4
#' alpha13 = 3.5
#' beta13 = 3
#' alpha23 = 2.5
#' beta23 = 1.5
#' tau = 0.3
#' p = 0.7
#' C_max = 10
#' alpha04 = 7
#' beta04 = 4
#' Result <- Estimation_unmasked(n, alpha12, beta12, alpha13, beta13, alpha23, beta23, tau, p, C_max, alpha04, beta04)
#' print(Result)
Estimation_unmasked<-function(n=200, alpha12 = 2, beta12 = 4, alpha13 = 3.5, beta13 = 3, alpha23 = 2.5, beta23 = 1.5, tau = 0.3, p = 0.7, C_max = 10, alpha04 = 7, beta04 = 4){


    Data <- data_generation_ver1(n=n, alpha12 = alpha12, beta12 = beta12, alpha13 = alpha13, beta13 = beta13, alpha23 = alpha23, beta23 = beta23, tau = tau, p = p, C_max = C_max, alpha04 = alpha04, beta04 = beta04)

    t1 <- Data$t1
    t2 <- Data$t2
    status12 = Data$status12
    status13 = Data$status13
    status23 = Data$status23
    n  = length(t1)

    ind = status12==1
    sum(ind)
    Cenrate_12 <- 1-sum(ind)/n

    ind = status13==1
    sum(ind)
    Cenrate_13 <- 1-sum(ind)/n

    ind = status23==1
    sum(ind)
    Cenrate_23 <- 1-sum(ind)/n

    # parameter of the Clayton copula
    phi <- 2/(1-tau)-2

    initial_parms <- c(alpha12, beta12, alpha13, beta13, alpha23, beta23, phi, p)


    Data_out <- nlm(log_like_func, initial_parms, hessian=TRUE)

    Data_out$n <- n
    Data_out$Cenrate_12 <- Cenrate_12
    Data_out$Cenrate_13 <- Cenrate_13
    Data_out$Cenrate_23 <- Cenrate_23

return(Data_out)
}






#' Simulation unmasked data
#'
#' @param B
#' @param n
#' @param alpha12
#' @param beta12
#' @param alpha13
#' @param beta13
#' @param alpha23
#' @param beta23
#' @param tau
#' @param p
#' @param C_max
#' @param alpha04
#' @param beta04
#' @return list of estimation outputs from 1:B
#' @export
#' @examples
#' B=3
#' n=200
#' alpha12 = 2
#' beta12 = 4
#' alpha13 = 3.5
#' beta13 = 3
#' alpha23 = 2.5
#' beta23 = 1.5
#' tau = 0.3
#' p = 0.7
#' C_max = 10
#' alpha04 = 7
#' beta04 = 4
#' Sim_result <- Simulation_unmasked(B, n, alpha12, beta12, alpha13, beta13, alpha23, beta23, tau, p, C_max, alpha04, beta04)
#' total_num <- length(Sim_result)
#'estimations <- c()
#'for (i in 1:total_num){
#'  estimation_temp <- Sim_result[[i]]$estimate
#'  estimations <- rbind(estimations, estimation_temp)
#'}
#'print(estimations)
Simulation_unmasked<-function(B=10, n=200, alpha12 = 2, beta12 = 4, alpha13 = 3.5, beta13 = 3, alpha23 = 2.5, beta23 = 1.5, tau = 0.3, p = 0.7, C_max = 10, alpha04 = 7, beta04 = 4){

  Data_out <-list()
  for (i in 1:B){
  Data_out_temp <- Estimation_unmasked(n, alpha12, beta12, alpha13, beta13, alpha23, beta23, tau, p, C_max, alpha04, beta04)
  Data_out[[i]] <- Data_out_temp
  }


  return(Data_out)
}
