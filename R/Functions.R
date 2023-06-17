#' Uniform Hazard Function
#'
#' @param a14 lower limit of the distribution.
#' @param b14 upper limit of the distribution.
#' @param t vector of quntiles.
#'
#' @return hazard function of uniform distribution.
#' @export
#'
#' @examples
#' hazardUnif(1,4,2)
hazardUnif = function(a14, b14, t){
  out = dunif(t, a14, b14)*(1-punif(t, a14, b14))
  return(out)
}



#' Cumulative Uniform Hazard
#'
#' @param a14 minimum value
#' @param b14 maximum value
#' @param t time to event t
#'
#' @return cumulative hazard function of uniform distribution
#' @export
#'
#' @examples
#' Cum_hazardUnif(1,4,2)
#' Cum_hazardUnif(1,4,0)  # outside of range gives 0
#' Cum_hazardUnif(1,4,5)  # outside of range gives 0
Cum_hazardUnif = function(a14, b14, t){
  if ((a14<t)&(t<b14)){
  out = log(dunif(t, a14, b14)*(1-punif(t, a14, b14)) )
  } else{
    out = 0 }
  return(out)
}








#' Nonparametric Estimated Cumulative Incidence Function
#'
#' @param time time
#' @param status12 indicator history of 1->2 transition
#' @param status13 indicator history of 1->3 transition
#' @param status14 indicator history of 1->4 transition
#'
#' @return Estimated Cumulative Incidence Function
#' @export
#'
#' @examples
#' library(dplyr)
#' library(ggplot2)
#' Data = data_generation()
#'
#'
#' status12 = Data$status12
#' status13 = Data$status12
#'
#' Data_new = Data %>% mutate(status14_new = status14+status04) %>% select(status14_new)
#' status14 = Data_new$status14_new
#' Data_results = CIF(Data$t1, status12, status13, status14)
#' ggplot(data=Data_results)+geom_point(aes(t, CIF1))


#' ### plot required...
CIF<-function(time,status12,status13,status14)
{
  K=3
  t=time

  S1<-status12
  S1[S1==1]<-1
  S2<-status13
  S2[S2==1]<-1
  S3<-status14
  S3[S3==1]<-1
  k=S1+S2+S3

  n=length(time)




  n.event.1<-status12

  n.event.2<-status13

  n.event.3<-status14
  n.event<-k


  data<-data.frame(cbind(t,n.event,n.event.1, n.event.2, n.event.3))

  data<-data[order(data$t),]
  n.risk<-rep(0,n)
  n.risk.1<-rep(0,n)
  n.risk.2<-rep(0,n)
  n.risk.3<-rep(0,n)


  n.risk[1]<-n
  n.risk.1[1]<-n
  n.risk.2[1]<-n
  n.risk.3[1]<-n

  n.event.1 = data$n.event.1
  n.event.2 = data$n.event.2
  n.event.3 = data$n.event.3



  n.risk=seq(1,n)
  n.risk=sort(n.risk,decreasing=TRUE)
  # n.risk =  cumsum(n.event.1+n.event.2+n.event.3)
  # n.risk = n.risk[200:1]

  CIF1<-rep(0,n)
  CIF2<-rep(0,n)
  CIF3<-rep(0,n)

  n.event = data$n.event

  surv<-rep(1,n)
  for (i in 1:n){
    pro<-1
    for (j in 1:i){
      pro<-(n.risk[j]-n.event[j])/(n.risk[j])
      surv[i]<-surv[i]*pro
    }
  }


  for (i in 1:n){
    SUM<-0
    for (j in 1:i){
      SUM<-surv[j]*(n.event.1[j]/n.risk[j])
      CIF1[i]<-CIF1[i]+SUM
    }
  }

  for (i in 1:n){
    SUM<-0
    for (j in 1:i){
      SUM<-surv[j]*(n.event.2[j]/n.risk[j])
      CIF2[i]<-CIF2[i]+SUM
    }
  }

  for (i in 1:n){
    SUM<-0
    for (j in 1:i){
      SUM<-surv[j]*(n.event.3[j]/n.risk[j])
      CIF3[i]<-CIF3[i]+SUM
    }
  }


  data_out<-data.frame(data, n.risk, surv, CIF1, CIF2, CIF3)

  return(data_out)
}









S10_fun = function(alpha12, beta12, alpha13, beta13, t){
  result = exp( -log(1+(t/alpha12)^beta12) -log(1+(t/alpha13)^beta13)   )
  return(result)
}


###
#' Title
#'
#' @param alphaj a scale parameter for transition 1->2.
#' @param betaj a shape parameter for transition 1->2.
#' @param alphak a scale parameter for transition 1->3.
#' @param betak a shape parameter for transition 1->3.
#' @param t time to event where t>0.
#'
#' @return
#' @export
#'
#' @examples f10j(1,2,1.3,2.3,5)
f10j<-function(alphaj, betaj, alphak, betak,t){
  ((betaj/alphaj)*(t/alphaj)^(betaj-1))/((1+(t/alphaj)^betaj)^2*(1+(t/alphak)^betak))
}




#' Title
#'
#' @param p
#'
#' @return
#' @export
#'
#' @examples
log.lik.fun<-function(p){


  ### M=1
  alpha12=p[1]
  beta12=p[2]


  #### M=2
  alpha13=p[3]
  beta2=p[4]


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
      f102<-f10j(alpha13, beta2, alpha12, beta12, t1[i])
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

