#' log-likelihood function without masked data
#'
#' @param p vector of parameters
#'
#' @return a value of log-likelihood
#' @export
#'
#' @examples
#' nlm(log_like_func, initial_value)
log_like_func<-function(p){


  ### M=1
  alpha12=p[1]
  beta12=p[2]


  #### M=2
  alpha13=p[3]
  beta13=p[4]


  ### M=3
  alpha23=p[5]
  beta23=p[6]



  ### parameter for Clayton Copula
  phi <- p[7]

  ### proportion
  prop <- p[8]

  ### M=2
  # psi2_cnot = p[9]




  answer=0
  for (i in 1:n){
    if (status12[i]==1){

      f101<-f10j(alpha12, beta12, alpha13, beta13, t1[i])
      F101 =  integrate(f10j, 0, t1[i], alphaj=alpha12, betaj=beta12, alphak=alpha13, betak=beta13)$value
      psi2 = integrate(f10j, 0, 20, alphaj=alpha12, betaj=beta12, alphak=alpha13, betak=beta13)$value

      s1<-1-F101/psi2

      #s1<-exp(-integrate(Hazard10,0, x1[i])$value)
      s3<-((1+(t1[i]/alpha12)^(beta12))^(-1)) * ((1+(t1[i]/alpha13)^(beta13))^(-1))
      ds1<--f101/psi2



      s2 <- exp(-(t2[i]/alpha23)^beta23)


      fun <- ((1-s1)**(-phi)-1)+((1-s2)**(-phi)-1)

      joints <- (fun+1)**(-1/phi)

      d1.joints <- (fun+1)**(-1/phi-1)*(1-s1)**(-phi-1)*(-ds1)

      if (status23[i]==1)
      {
        ## maybe check it here
        ds2 <- -beta23*(t2[i]^(beta23-1)/alpha23^beta23)*exp(-(t2[i]/alpha23)^beta23)
        d2.joints <- (fun+1)**(-1/phi-1)*(1-s2)**(-phi-1)*(-ds2)

        dd.joints <- d1.joints*d2.joints*(fun+1)**(1/phi)*(1+phi)


      }

    }else if(status13[i]==1){
      f102<-f10j(alpha13, beta13, alpha12, beta12, t1[i])
    }else if((status12[i]==0)&(status13[i]==0)){

      S10<-((1+(t1[i]/alpha12)^(beta12))^(-1)) * ((1+(t1[i]/alpha13)^(beta13))^(-1))


    }
    if ((status12[i]==1)&(status23[i]==1)){answer <- answer-log(psi2)-log(prop)-log(dd.joints)}
    if ((status12[i]==1)&(status23[i]==0)) {answer <- answer-log(prop* psi2*(-ds1)-prop* psi2*d1.joints)}
    if (status13[i]==1){answer=answer-log(prop)-log(f102)}
    if ((status12[i]==0)&(status13[i]==0)){answer=answer-log(prop*S10+(1-prop))}

  }
  return(answer)



}



# likelihood function for masked causes
log.lik.fun.em.adapIntegrate<-function(p){


  ### M=1
  alpha1=p[1]
  beta1=p[2]


  #### M=2
  alpha2=p[3]
  beta2=p[4]

  ### M=3
  lambda23=p[5]
  k23=p[6]

  # ### M=4
  # lambda04=p[7]
  # k04=p[8]

  ### parameter for Clayton Copula
  phi <- p[7]

  ### proportion
  prop <- p[8]

  ###
  P_gm_gb = p[9]
  P_gm_gc = p[10]



  answer=0
  for (i in 1:n){
    if (status12[i]==1){
      # pie1 = integrate(f10j, 0, 20, alphaj=alpha1, betaj=beta1, alphak=alpha2, betak=beta2)$value
      pie1 = adaptIntegrate(f10j, 0, 20, alphaj=alpha1, betaj=beta1, alphak=alpha2, betak=beta2)$integral

      f101<-f10j(alpha1, beta1, alpha2, beta2, t1[i])
      # F101 =  integrate(f10j, 0, t1[i], alphaj=alpha1, betaj=beta1, alphak=alpha2, betak=beta2)$value
      F101 =  adaptIntegrate(f10j, 0, t1[i], alphaj=alpha1, betaj=beta1, alphak=alpha2, betak=beta2)$integral


      s1<-1-F101/pie1
      #s1<-exp(-integrate(Hazard10,0, t1[i])$value)
      s3<-((1+(t1[i]/alpha1)^(beta1))^(-1)) * ((1+(t1[i]/alpha2)^(beta2))^(-1))
      ds1<--f101/pie1


      s2 <- exp(-(t2[i]/lambda23)^k23)


      fun <- ((1-s1)**(-phi)-1)+((1-s2)**(-phi)-1)

      joints <- (fun+1)**(-1/phi)

      d1.joints <- (fun+1)**(-1/phi-1)*(1-s1)**(-phi-1)*(-ds1)

      if (status23[i]==1)
      {
        ds2 <- -k23*(t2[i]^(k23-1)/lambda23^k23)*exp(-(t2[i]/lambda23)^k23)
        d2.joints <- (fun+1)**(-1/phi-1)*(1-s2)**(-phi-1)*(-ds2)

        dd.joints <- d1.joints*d2.joints*(fun+1)**(1/phi)*(1+phi)


      }

    }else if((status13[i]==1)&(gamma_gm[i]==1)){
      # pie1 = adaptIntegrate(f10j, lowerLimit=0, upperLimit=max(t1), alphaj=alpha1, betaj=beta1, alphak=alpha2, betak=beta2)$integral

      f102star<-f10j(alpha2, beta2, alpha1, beta1, t1[i])
    }else if((status13[i]==1)&(gamma_gm[i]==0)){
      # pie1 = adaptIntegrate(f10j, lowerLimit=0, upperLimit=max(t1), alphaj=alpha1, betaj=beta1, alphak=alpha2, betak=beta2)$integral

      f102star<-f10j(alpha2, beta2, alpha1, beta1, t1[i])
    }else if((status13[i]==-1)&(gamma_gm[i]==1)){
      # pie1 = adaptIntegrate(f10j, lowerLimit=0, upperLimit=max(t1), alphaj=alpha1, betaj=beta1, alphak=alpha2, betak=beta2)$integral

      f102star<-f10j(alpha2, beta2, alpha1, beta1, t1[i])
      S10<-((1+(t1[i]/alpha1)^(beta1))^(-1)) * ((1+(t1[i]/alpha2)^(beta2))^(-1))

    }else if((status12[i]==0)&(status13[i]==0)&(gamma_gm[i]==0)){

      S10<-((1+(t1[i]/alpha1)^(beta1))^(-1)) * ((1+(t1[i]/alpha2)^(beta2))^(-1))

    }else if((status12[i]==0)&(status13[i]==0)&(gamma_gm[i]==1)){

      S10<-((1+(t1[i]/alpha1)^(beta1))^(-1)) * ((1+(t1[i]/alpha2)^(beta2))^(-1))

    }else if((status12[i]==0)&(status13[i]==-1)&(gamma_gm[i]==1)){

      S10<-((1+(t1[i]/alpha1)^(beta1))^(-1)) * ((1+(t1[i]/alpha2)^(beta2))^(-1))

    }

    if ((status12[i]==1)&(status23[i]==1)){answer <- answer-log(prop)-log(pie1)-log(dd.joints)}
    if ((status12[i]==1)&(status23[i]==0)) {answer <- answer-log(prop*pie1*(-ds1)-prop*pie1*d1.joints)}
    if ((status13[i]==1)&(gamma_gm[i]==1)){answer=answer-log(prop*f102star)-log(P_gm_gb)}
    if ((status13[i]==-1)&(gamma_gm[i]==1)){answer=answer-Pi_gb_gm_hat[i]*log(prop*f102star)-Pi_gb_gm_hat[i]*log(P_gm_gb)}
    if ((status13[i]==1)&(gamma_gm[i]==0)){answer=answer-log(prop*f102star)-log(1-P_gm_gb)}

    if ((status12[i]==0)&(status13[i]==0)&(gamma_gm[i]==0)&(status_d[i]==1)){answer=answer-log(prop*S10+(1-prop))-log(1-P_gm_gc)}
    if ((status12[i]==0)&(status13[i]==0)&(gamma_gm[i]==1)&(status_d[i]==1)){answer=answer-log(prop*S10+(1-prop))-log(P_gm_gc+0.001)}
    if ((status12[i]==0)&(status13[i]==-1)&(gamma_gm[i]==1)){answer=answer-(1-Pi_gb_gm_hat[i])*log(prop*S10+(1-prop))-(1-Pi_gb_gm_hat[i])*log(P_gm_gc+0.001)}

  }
  return(answer)
}







#' Constructor function that creates a negative log-likelihood function with unmaksed data
#'
#' @param data generated data
#' @param fixed parameters
#'
#' @return negative log-likelihood function
#' @export
#'
#' @examples
#' Data <- data_generation_ver1(n=n, alpha12 = alpha12, beta12 = beta12, alpha13 = alpha13, beta13 = beta13, alpha23 = alpha23, beta23 = beta23, tau = tau, p = p, C_max = C_max, alpha04 = alpha04, beta04 = beta04)
#' nLL <- make.NegLogLik(Data)
#' phi <- 2/(1-tau)-2 # parameter of the Clayton copula
#' initial_parms <- c(alpha12, beta12, alpha13, beta13, alpha23, beta23, phi, p)
#' estimates <- nlm(nLL, initial_parms)$estimate
#' print(estimates)
make.NegLogLik<-function(data, fixed = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)){

  params <- fixed

  function(p){
    ### M=1
    alpha12=p[1]
    beta12=p[2]


    #### M=2
    alpha13=p[3]
    beta13=p[4]


    ### M=3
    alpha23=p[5]
    beta23=p[6]



    ### parameter for Clayton Copula
    phi <- p[7]

    ### proportion
    prop <- p[8]

    ### M=2
    # psi2_cnot = p[9]

    n <- length(data$t1)
    t1 <- data$t1
    t2 <- data$t2
    status12 = data$status12
    status13 = data$status13
    status23 = data$status23
    n  = length(t1)

    max_t1 <- max(t1)

    answer=0
    for (i in 1:n){
      if (status12[i]==1){

        f101<-f10j(alpha12, beta12, alpha13, beta13, t1[i])



        F101 <- cubature::adaptIntegrate(f10j, 0, t1[i], alphajk=alpha12, betajk=beta12, alphajl=alpha13, betajl=beta13)$integral
        psi2 <- cubature::adaptIntegrate(f10j, 0, max_t1, alphajk=alpha12, betajk=beta12, alphajl=alpha13, betajl=beta13)$integral

        s1 <- 1-F101/psi2

        #s1<-exp(-integrate(Hazard10,0, x1[i])$value)
        s3<-((1+(t1[i]/alpha12)^(beta12))^(-1)) * ((1+(t1[i]/alpha13)^(beta13))^(-1))
        ds1<--f101/psi2


        s2 <- exp(-(t2[i]/alpha23)^beta23)

        fun <- ((1-s1)**(-phi)-1)+((1-s2)**(-phi)-1)

        joints <- (fun+1)**(-1/phi)

        d1.joints <- (fun+1)**(-1/phi-1)*(1-s1)**(-phi-1)*(-ds1)

        if (status23[i]==1)
        {
          ## maybe check it here
          ds2 <- -beta23*(t2[i]^(beta23-1)/alpha23^beta23)*exp(-(t2[i]/alpha23)^beta23)
          d2.joints <- (fun+1)**(-1/phi-1)*(1-s2)**(-phi-1)*(-ds2)

          dd.joints <- d1.joints*d2.joints*(fun+1)**(1/phi)*(1+phi)


        }

      }else if(status13[i]==1){
        f102<-f10j(alpha13, beta13, alpha12, beta12, t1[i])
      }else if((status12[i]==0)&(status13[i]==0)){

        S10<-((1+(t1[i]/alpha12)^(beta12))^(-1)) * ((1+(t1[i]/alpha13)^(beta13))^(-1))


      }
      if ((status12[i]==1)&(status23[i]==1)){answer <- answer-log(psi2)-log(prop)-log(dd.joints)}
      if ((status12[i]==1)&(status23[i]==0)) {answer <- answer-log(prop* psi2*(-ds1)-prop* psi2*d1.joints)}
      if (status13[i]==1){answer=answer-log(prop)-log(f102)}
      if ((status12[i]==0)&(status13[i]==0)){answer=answer-log(prop*S10+(1-prop))}

    }
    return(answer)


  }
}




#' Constructor function that creates a negative log-likelihood function with unmaksed data first stage
#'
#' @param data generated data
#' @param fixed parameters
#'
#' @return negative log-likelihood function for first stage
#' @export
#'
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
#' Data <- data_generation_ver1(n=n, alpha12 = alpha12, beta12 = beta12, alpha13 = alpha13, beta13 = beta13, alpha23 = alpha23, beta23 = beta23, tau = tau, p = p, C_max = C_max, alpha04 = alpha04, beta04 = beta04)
#' nLL <- make.NegLogLik_first_stage(Data)
#' initial_parms <- c(alpha12, beta12, alpha13, beta13, p)
#' nlm(nLL, initial_parms)$estimate
make.NegLogLik_first_stage<-function(data, fixed = c(FALSE, FALSE, FALSE, FALSE, FALSE)){

  params <- fixed

  function(p){
    ### M=1
    alpha12=p[1]
    beta12=p[2]


    #### M=2
    alpha13=p[3]
    beta13=p[4]


    ### proportion
    prop <- p[5]



    n <- length(data$t1)
    t1 <- data$t1
    t2 <- data$t2
    status12 = data$status12
    status13 = data$status13
    status23 = data$status23
    n  = length(t1)

    max_t1 <- max(t1)

    answer=0
    for (i in 1:n){
      if (status12[i]==1){

        f101<-f10j(alpha12, beta12, alpha13, beta13, t1[i])


      }else if(status13[i]==1){
        f102<-f10j(alpha13, beta13, alpha12, beta12, t1[i])
      }else if((status12[i]==0)&(status13[i]==0)){

        S10<-((1+(t1[i]/alpha12)^(beta12))^(-1)) * ((1+(t1[i]/alpha13)^(beta13))^(-1))


      }
      if ((status12[i]==1)){answer <- answer-log(prop)-log(f101)}
      if (status13[i]==1){answer=answer-log(prop)-log(f102)}
      if ((status12[i]==0)&(status13[i]==0)){answer=answer-log(prop*S10+(1-prop))}

    }
    return(answer)

  }
}



#' Constructor function that creates a negative log-likelihood function with unmaksed data second stage
#'
#' @param data generated data
#' @param fixed parameters
#'
#' @return negative log-likelihood function for second stage
#' @export
#'
#' @examples
#' Data <- data_generation_ver1(n=n, alpha12 = alpha12, beta12 = beta12, alpha13 = alpha13, beta13 = beta13, alpha23 = alpha23, beta23 = beta23, tau = tau, p = p, C_max = C_max, alpha04 = alpha04, beta04 = beta04)
#' cumF101 <- setRefClass("cumulative",
#' fields = list(value = "numeric", alpha12 = "numeric", beta12 = "numeric", alpha13 = "numeric", beta13 = "numeric"),
#' methods = list(
#'   integrate = function(x) {
#'     value <<- cubature::adaptIntegrate(Vectorize(f10j), lowerLimit=0, upperLimit=x, alphajk=alpha12, betajk=beta12, alphajl=alpha13, betajl=beta13)$integral
#'     return(value)
#'   }
#' )
#' )
#' nLL <- make.NegLogLik_first_stage(Data)
#' initial_parms <- c(alpha12, beta12, alpha13, beta13, p)
#' estimates <- nlm(nLL, initial_parms)$estimate
#' alpha12hat = estimates[1]
#' beta12hat = estimates[2]
#' alpha13hat = estimates[3]
#' beta13hat = estimates[4]
#' cumF101_hat <- cumF101$new(alpha12 = alpha12hat, beta12 = beta12hat, alpha13 = alpha13hat, beta13=beta13hat)
#' t1 <- Data$t1
#' F101_hat <- sapply(t1, cumF101_hat$integrate)
#' F101.star.hat.new <- F101_hat
#' pi1<-hazardLL(alpha12hat,beta12hat,t1)/(hazardLL(alpha12hat,beta12hat,t1)+hazardLL(alpha13hat,beta13hat,t1))
#' for (i in 1:length(F101_hat)){
#'F101.star.hat.new[i]=F101_hat[i]/pi1[i]
#'}
#' Data$F101.star.hat <-F101.star.hat.new
#' nLL2 <- make.NegLogLik_second_stage(Data)
#' phi <- 2/(1-tau)-2 # parameter of the Clayton copula
#' initial_parms <- c(0.7,1.3,0.8)
#' estimates2<-nlm(nLL2, initial_parms)$estimate
#' print(estimates)
#' print(estimates2)
make.NegLogLik_second_stage<-function(data, fixed = c(FALSE, FALSE, FALSE)){

  params <- fixed

  function(p){
    # 2 -> 3
    # params[!fixed] <- p
    alpha23 <- params[1]
    beta23 <- params[2]
    phi <- params[3]

    n <- length(data$t1)
    t1 <- data$t1
    t2 <- data$t2
    status12 = data$status12
    status13 = data$status13
    status23 = data$status23
    n  = length(t1)
    F101.star.hat.new <- data$F101.star.hat

    answer <- 0
    for (i in 1:n)
    {

      if (status12[i]==1)
      {

        s2 <- exp(-(t2[i]/alpha23)^beta23)
        f2 <-1-s2

        fun <- ((F101.star.hat.new[i])**(-phi)-1)+((1-s2)**(-phi)-1)

        joints <- (fun+1)**(-1/phi)   ### Copula

        d1.joints <- (fun+1)**(-1/phi-1)*(F101.star.hat.new[i])**(-phi-1)

        if (status23[i]==1)
        {
          ds2 <- -beta23*(t2[i]^(beta23-1)/alpha23^beta23)*exp(-(t2[i]/alpha23)^beta23)

          dd.joints <- F101.star.hat.new[i]^(-phi-1)*f2^(-phi-1)*(fun+1)**(-1/phi-2)*(1+phi)*(ds2)

        }
      }



      if ((status12[i]==1)&(status23[i]==1)) {answer <- answer-log(dd.joints+0.001)}
      if ((status12[i]==1)&(status23[i]==0)) {answer <- answer-log(1-d1.joints)}



    }

    return(answer)

  }
}




#' Title
#'
#' @param p parameters of 2->3 and clayton copula
#'
#' @return negative loglikehood values
#' @export
#'
#' @examples
#'
#' f101<-sapply(t1, f10j, alphajk=alpha12hat, betajk=beta12hat, alphajl=alpha13hat, betajl=beta13hat)
#f101<-f10j(alpha1hat, beta1hat, alpha2hat, beta2hat, x1[i])
#' cumF101<-function(t)
#'{
#'  return=cubature::adaptIntegrate(f10j, lowerLimit=0, upperLimit=t, alphajk=alpha12hat, betajk=beta12hat, alphajl=alpha13hat, betajl=beta13hat)$integral
#'}
#'
#' F101<-sapply(t1,cumF101)
#'F101.star.hat.new=F101
#'pi1<-rep(0,length(F101))
#'
#'for (i in 1:n){
  #'  pi1[i]<-cubature::adaptIntegrate(f10j, lowerLimit=0, upperLimit=t1[i], alphajk=alpha12hat, betajk=beta12hat, alphajl=alpha13hat, betajl=beta13hat)$integral
  #'}
#'
#' for (i in 1:length(F101)){
#'  F101.star.hat.new[i]=F101[i]/pi1[i]
#' }
#' status12 <- Data$status12
#' status23 <- Data$status12
#' F101.star.hat.new <- Data$F101.star.hat
#' t2 <- Data$t2
#' results2 <- nlm(log.lik.clayton.np,c(0.7,1.3,0.8))
#' print(results2)
log_lik_second_stage <- function(p)
{

  ### M=3
  alpha23=p[1]
  beta23=p[2]



  ### parameter for Clayton Copula
  phi <- p[3]


  n=length(status12)
  answer <- 0
  for (i in 1:n)
  {

    if (status12[i]==1)
    {






      s2 <- exp(-(t2[i]/alpha23)^beta23)
      f2 <-1-s2

      fun <- ((F101.star.hat.new[i])**(-phi)-1)+((1-s2)**(-phi)-1)

      joints <- (fun+1)**(-1/phi)   ### Copula

      d1.joints <- (fun+1)**(-1/phi-1)*(F101.star.hat.new[i])**(-phi-1)

      if (status23[i]==1)
      {
        ds2 <- beta23*(t2[i]^(beta23-1)/alpha23^beta23)*exp(-(t2[i]/alpha23)^beta23)

        dd.joints <- F101.star.hat.new[i]^(-phi-1)*f2^(-phi-1)*(fun+1)**(-1/phi-2)*(1+phi)*(ds2)


      }
    }



    if ((status12[i]==1)&(status23[i]==1)) {answer <- answer-log(dd.joints)}
    if ((status12[i]==1)&(status23[i]==0)) {answer <- answer-log(1-d1.joints)}



  }
  answer
}



