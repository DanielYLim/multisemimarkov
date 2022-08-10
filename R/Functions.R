
#' Log-Logistic hazard function
#'
#' @param alpha1 a scale parameter.
#' @param beta1  a shape parameter.
#' @param t time to event where t>0.
#'
#' @return hazard function of t with parameters alpha1 adn beta1.
#' @export
#'
#' @examples
#' hazardLL(1,2,1.3)
#' \dontrun{
#' hazardLL(1,2,2.3)
#' }
hazardLL<-function(alpha1,beta1,t){
	((beta1/alpha1)*(t/alpha1)^(beta1-1))/(1+(t/alpha1)^beta1)
}



#' Weibull hazard function
#'
#' @param k a shape parameter.
#' @param lambda a scale parameter.
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
hazardWei = function(k, lambda, t){
	(k/lambda)*(t/lambda)**(k-1)
}


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
  out = ((t-a14)*(b14-a14)-(t-a14)**2)/( b14-a14 )**2
  return(out)
}










#' Data Generation Function for a Illness-Death model with Masked Causes
#'
#' @param n number of observations.
#' @param alpha12 log-logstic scale parameter; 1->2 transition.
#' @param beta12 log-logstic shape parameter; 1->2 transition.
#' @param alpha13 log-logstic scale parameter; 1->3 transition.
#' @param beta13 log-logstic shape parameter; 1->3 transition.
#' @param lambda23 Weibull scale parameter; 2->3 transition.
#' @param k23 Weibull shape parameter; 2->3 transition.
#' @param lambda04 Weibull scale parameter; 0->4 transition.
#' @param k04 Weibull shape parameter; 0->4 transition.
#' @param a14 Uniform lower limit; 1->4 transition with non cured proportion.
#' @param b14 Uniform upper limit; 1->4 transition with non cured proportion.
#' @param tau Pre-specified upper threshold value for cured individual.
#' @param p Probability of not-being cured.
#' @param Pg_mb Probability of masked given that it was 1->3 transition.
#' @param Pg_mc Probability of masked given that it was 1->4 transition.
#' @param C_max Maximum value of censoring time.
#' @param tau_0 Pre-specified lower threshold value for not being cured individual.
#'
#' @return Data frame
#' @export
#'
#' @examples
#' require(dplyr)
#' Data_temp = data_generation(n=200, alpha12=2.5, beta12=2, alpha13=5, beta13=2.5, lambda23=1, k23=2, lambda04=10, k04=4, a14=4, b14=2, tau=0.3, p=0.7, Pg_mb=0.3, Pg_mc=0.2, C_max=10, tau_0=2)
#' head(Data_temp)
#'
#'

data_generation = function(n=200, alpha12=2.5, beta12=2, alpha13=5, beta13=2.5, lambda23=1, k23=2, lambda04=4, k04=2, a14=1, b14=4, tau=0.3, p=0.7, Pg_mb=0.3, Pg_mc=0.2, C_max=10, tau_0=2){


  gett<-function(t){
    Sol<-log(1-u)+log(1+(t/alpha12)^(beta12))+log(1+(t/alpha13)^(beta13))+log( ((t-a14)*(b14-a14)-(t-a14)**2)/( b14-a14 )**2   )
    return(Sol)
  }

  phi <- 2/(1-tau)-2

  status12<-rep(1,n)
  status13<-rep(1,n)
  status23<-rep(1,n)
  status14<-rep(1,n)
  status04<-rep(1,n)

  x12<-rep(0,n)
  x13<-rep(0,n)
  x23<-rep(1,n)
  x14<-rep(0,n)
  x04<-rep(0,n)


  W<-rbinom(n,1,p)


  ### Censoring times
  C<-runif(n,0,C_max)



  u1<-runif(n,0,1)
  u4<-runif(n,0,1)


  for (i in 1:n){
    if (W[i]==1){
      u<-u1[i]
      t<-uniroot(gett,c(0,C_max+15),extendInt="yes")$root
      pi1<-hazardLL(alpha12,beta12,t)/(hazardLL(alpha12,beta12,t)+hazardLL(alpha13,beta13,t)+hazardUnif(a14, b14, t)) 	 ## prob of having t with 1->2
      pi2<-hazardLL(alpha13,beta13,t)/(hazardLL(alpha12,beta12,t)+hazardLL(alpha13,beta13,t)+hazardUnif(a14, b14, t))  	## prob of having t with 1->3
      pi3 = 1-(pi1+pi2)
      prob = c(pi1, pi2, pi3)

      S<-rmultinom(1,1,prob)

      if (S[1,1]==1){
        x12[i]=t
        x13[i]=0
        x14[i]=0
        x04[i]=0
        status12[i]=1
        status13[i]=0
        status14[i]=0
        status04[i]=0
      }else if(S[2,1]==1){
        x12[i]=0
        x13[i]=t
        x14[i]=0
        x04[i]=0
        status12[i]=0
        status13[i]=1
        status14[i]=0
        status04[i]=0
      }else if(S[3,1]==1){
        x12[i]=0
        x13[i]=0
        x14[i]=t
        x04[i]=0
        status12[i]=0
        status13[i]=0
        status14[i]=1
        status04[i]=0
      }



    }else if (W[i]==0){
      x12[i]=0
      x13[i]=0
      x14[i]=0
      x04[i]=lambda04*(-log(1-u1[i]))^(1/k04)
      status12[i]=0
      status13[i]=0
      status14[i]=0
      status04[i]=1

    }
  }



  # generate data from the Clayton copula with given marginal distributions
  for (i in 1:n)
  {
    if (status12[i]==1)
    {

      u2 <- runif(1,min=0,max=1)
      #u2 = 1-pweibull(rweibull(1,k23,lambda23), k23, lambda23)
      #x23[i] <- sqrt(-log(1-(1-u1[i]**(-phi)*(1-u2**(-phi/(1+phi))))**(-1/phi)))   ##### scale parameter 1, shape parameter 2
      x23[i] = lambda23*(-log(1-(1-u1[i]**(-phi)*(1-u2**(-phi/(1+phi))))**(-1/phi)))**(1/k23)
    }else if (status12[i]==0)
    {
      status23[i]=0
      x23[i]=0
    }
  }



  for (i in 1:n){
    if ((x12[i]>C[i])&(status12[i]==1)){
      x12[i]=C[i]
      x13[i]=0
      x14[i]=0
      x23[i]=0
      status12[i]=0
      status23[i]=0
      status13[i]=0
      status14[i]=0
    }
    if ((status12[i]==1)&(x12[i]+x23[i]>C[i]))
    {
      status23[i] <- 0
      x23[i] <- C[i]-x12[i]
    }
    if ((x13[i]>C[i])&(status13[i]==1)){
      x13[i]=C[i]
      status13[i]=0
    }
    if ((x14[i]>C[i])&(status14[i]==1)){
      x14[i]=C[i]
      status14[i]=0
    }
    if ((x04[i]>C[i])&(status04[i]==1)){
      x04[i]=C[i]
      status04[i]=0
    }
  }




  t1<-x12+x13+x14+x04
  t2<-x23


  # status_c = status12+status13
  # status_c = (status_c==0)

  gamma_gm<-rep(0,n)


  ### Missing Data Generation
  status_a<- status12
  status_b<- status13
  status_c<- status14+status04


  Index<-which(status_c==1)

  temp_num1=0
  temp_num2=0
  for (i in 1:length(Index)){
    NUM<-rbinom(1,1,Pg_mc)
    if (NUM==1){
      status_b[Index[i]]=-1
      status_c[Index[i]]=-1
      temp_num1=temp_num1+1
    }
  }

  Index<-which(status13 == 1)

  for (i in 1:length(Index)){
    NUM<-rbinom(1,1,Pg_mb)
    if (NUM==1){
      status_b[Index[i]]=-1
      status_c[Index[i]]=-1
      temp_num2=temp_num2+1
    }
  }

  NUM<-which(status_b == -1)
  gamma_gm[NUM]=1

  NUM<-which(status_c == -1)
  gamma_gm[NUM]=1





  #### Data imputation procedure
  #### about 50 % come back (e.g.)

  ########## Missing Data is generated

  status13new<-status_b
  status13new[status13new==-1]<-0
  status1new<-status12+status13new

  tau_1=max(t1[which(status1new==1)]) ### point

  # tau_0 = 2





  ## Assign the missing indicator to the cured to death comparing with the cutoff point
  temp4 = 0
  ### Imputation 1
  Number<-which(status_c==-1)
  Number2<-t1[Number]>tau_1
  Number3<-which(Number2==TRUE)
  Number4<-Number[Number3]
  temp4=length(Number4)

  status_c[Number4]<-1
  status_b[Number4]<-0


  ### Imputation 2
  Number<-which(status_b==-1)
  Number2<-t1[Number]<tau_0
  Number3<-which(Number2==TRUE)
  Number4<-Number[Number3]
  temp5=length(Number4)

  status_b[Number4]<-1
  status_c[Number4]<-0


  Data = data.frame(t1, t2, C, status_a, status_b, status_c, status12, status23, status13, status14, status04, gamma_gm)

  ### Order the data by t1
  Data <-Data[order(Data$t1),]
  return(Data)
}



S10_fun = function(alpha12, beta12, alpha13, beta13, t){
  result = exp( -log(1+(t/alpha12)^beta12) -log(1+(t/alpha13)^beta13)   )
  return(result)
}


###
f10j<-function(alphaj, betaj, alphak, betak,t){
  ((betaj/alphaj)*(t/alphaj)^(betaj-1))/((1+(t/alphaj)^betaj)^2*(1+(t/alphak)^betak))
}


f101_fun<-function(alpha12, beta12, alpha13, beta13, t){
  h101 = hazardLL(alpha12, beta12,t)
  S10 = S10_fun(alpha12, beta12, alpha13, beta13, t)
  return(h101 * S10 )
}

# f101_fun<-function(alpha12, beta12, alpha13, beta13, lambda_c, k_c, t){
#   ((beta12/alpha12)*(t/alpha12)^(beta12-1))/(1+(t/alpha12)^beta12) * exp( -log(1+(t/alpha12)^beta12) -log(1+(t/alpha13)^beta13)  - (t/lambda_c)^(k_c)  )
# }


f102_fun<-function(alpha12, beta12, alpha13, beta13, t){
  h102 = hazardLL(alpha13, beta13,t)
  S10 = S10_fun(alpha12, beta12, alpha13, beta13, t)
  return(h102 * S10 )
}


