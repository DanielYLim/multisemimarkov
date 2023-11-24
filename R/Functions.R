#' Uniform Hazard Function
#'
#' @param a14 lower limit of the distribution.
#' @param b14 upper limit of the distribution.
#' @param t vector of quantiles.
#'
#' @return a value of hazard using uniform distribution.
#' @export
#'
#' @examples
#' hazardUnif(1,4,2)
hazardUnif = function(a14, b14, t){
  if ((a14<t)&(t<b14)){
    out <- dunif(t, a14, b14)/(1-punif(t, a14, b14))
  }else{
    out <- 0
  }

  return(out)
}



#' Cumulative Uniform Hazard
#'
#' @param a14 minimum value
#' @param b14 maximum value
#' @param t time to event t
#'
#' @return a value of cumulative hazard using uniform distribution
#' @export
#'
#' @examples
#' Cum_hazardUnif(1,4,2)
#' Cum_hazardUnif(1,4,0)  # outside of range gives 0
#' Cum_hazardUnif(1,4,5)  # outside of range gives 0
Cum_hazardUnif = function(a14, b14, t){
  if ((a14<t)&(t<b14)){
  out = -log((1-punif(t, a14, b14)) )
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
#' @return A data frame contains estimated cumulative incidence function
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









S10_fun = function(alphajk, betajk, alphajl, beta13, t){
  result = exp( -log(1+(t/alphajk)^betajk) -log(1+(t/alphajl)^beta13)   )
  return(result)
}


###
#' Title
#'
#' @param alphajk a scale parameter for transition j->k.
#' @param betajk a shape parameter for transition j->2k.
#' @param alphajl a scale parameter for transition j->l.
#' @param betajl a shape parameter for transition j->l.
#' @param t time to event where t>0.
#'
#' @return a value of hazard
#' @export
#'
#' @examples
#' f10j(1,2,1.3,2.3,5)
f10j<-function(alphajk, betajk, alphajl, betajl,t){
  ((betajk/alphajk)*(t/alphajk)^(betajk-1))/((1+(t/alphajk)^betajk)^2*(1+(t/alphajl)^betajl))
}




cumF101 <- setRefClass("cumulative",
                       fields = list(value = "numeric", alpha12 = "numeric", beta12 = "numeric", alpha13 = "numeric", beta13 = "numeric"),
                       methods = list(
                          integrate = function(x) {
                           value <<- cubature::adaptIntegrate(Vectorize(f10j), lowerLimit=0, upperLimit=x, alphajk=alpha12, betajk=beta12, alphajl=alpha13, betajl=beta13)$integral
                           return(value)
                         }
                       )
)




