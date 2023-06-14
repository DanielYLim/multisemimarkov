

#' Title
#'
#' @param n number of observations.
#' @param alpha12 log-logstic scale parameter; 1->2 transition.
#' @param beta12 log-logstic shape parameter; 1->2 transition.
#' @param alpha13 log-logstic scale parameter; 1->3 transition.
#' @param beta13 log-logstic shape parameter; 1->3 transition.
#' @param lambda23 Weibull scale parameter; 2->3 transition.
#' @param k23 Weibull shape parameter; 2->3 transition.
#' @param tau Kendall's tau; dependency measure
#' @param p Probability of not-being cured.
#' @param C_max Maximum value of censoring time.
#'
#' @return
#' @export
#'
#' @examples
#' require(dplyr)
#' Data_temp = data_generation_ver1(n=1000)
#' head(Data_temp)
#' Data_temp %>% select(t1, t2, status12, status23, status13) %>% filter(status23==1)
data_generation_ver1 = function(n, alpha12=1, beta12=1, alpha13=2, beta13=2, lambda23=1, k23=2, tau=0.2, p=0.7, C_max=10){

  status12<-rep(1,n)
  status23<-rep(1,n)
  status13<-rep(1,n)

  x12<-rep(0,n)
  x23<-rep(0,n)
  x13<-rep(0,n)


  W<-rbinom(n,1,p)

  phi <- 2/(1-tau)-2
  ### Censoring times

  C<-runif(n,0,10)
  u1<-runif(n,0,1)
  u4<-runif(n,0,1)

  gett<-function(t){
    Sol<-log(1-u)+log(1+(t/alpha12)^(beta12))+log(1+(t/alpha13)^(beta13))
    return(Sol)
  }

  for (i in 1:n){
    if (W[i]==1){
      u<-u1[i]
      t<-uniroot(gett,c(0,20),extendInt="yes")$root
      pi1<-hazardLL(alpha12,beta12,t)/(hazardLL(alpha12,beta12,t)+hazardLL(alpha13,beta13,t))
      S<-rbinom(1,1,pi1)
      if (S==1){
        x12[i]=t
        x13[i]=0
        status12[i]=1
        status13[i]=0
      }else if(S==0){
        x12[i]=0
        x13[i]=t
        status12[i]=0
        status13[i]=1
      }
    }else if (W[i]==0){
      x12[i]=0
      x13[i]=0
      status12[i]=0
      status23[i]=0
      status13[i]=0
    }
  }



  # generate data from the Clayton copula with given marginal distributions
  for (i in 1:n)
  {
    if (status12[i]==1)
    {

      u3 <- runif(1,min=0,max=1)
      x23[i] <- sqrt(-log(1-(1-u1[i]**(-phi)*(1-u3**(-phi/(1+phi))))**(-1/phi)))

    }else if (status12[i]==0)
    {
      status23[i]=0
      x23[i]=0
    }
  }





  for (i in 1:n){
    if ((x12[i]>C[i])&(status12[i]==1)){
      x12[i]=C[i]
      x23[i]=0
      status12[i]=0
      status23[i]=0
    }
    if ((x23[i]>C[i])&(status23[i]==1)){
      x23[i]=C[i]
      status23[i]=0
    }
    if ((status12[i]==1)&(x12[i]+x23[i]>C[i]))
    {
      status23[i] <- 0
      x23[i] <- C[i]-x12[i]
    }
  }
  t1<-x12+x13
  t2<-x23

  ind = t1==0
  t1[ind] = C[ind]


  Data = cbind(t1,t2,status12, status23, status13)
  Data = as.data.frame(Data)
  return(Data)
}
