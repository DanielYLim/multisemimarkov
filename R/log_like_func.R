#' Title
#'
#' @param p
#'
#' @return
#' @export
#'
#' @examples
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


