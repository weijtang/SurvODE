args = commandArgs(trailingOnly=TRUE) # sample size

########################## covariance matrix 
AR1 <- function(tau, m) {
  if(m==1) {R <- 1}
  if(m > 1) {
    R <- diag(1, m)
    for(i in 1:(m-1)) {
      for(j in (i+1):m) {
        R[i,j] <- R[j,i] <- tau^(abs(i-j))
      }
    }
  }
  return(R)
}
#######################################


library(survival)
library(tmvtnorm)
library(splines2)
p=5
N=as.integer(args[1])

loop=0
pb = txtProgressBar(style = 3)
nloop=1000
while (loop<nloop) {
  loop=loop+1
  set.seed(loop)
  ############generate data########################
  Sigma_z1<-AR1(0.6,p)
  
  z= rtmvnorm(N, mean=rep(0,p), sigma=Sigma_z1,
              lower=rep(-2, length = p),
              upper=rep( 2, length = p))
  
  U=runif(N, 0,1)
  
  pre_time=rep(0, N)
  for (i in 1:(N)) {
    f=function(t) {
      integrand <- function(x) {0.5*exp(z[i,1]-z[i,3]+sin(3*pi*x/4)*(x<3)*z[i,2]-z[i,4]+z[i,5])}
      
      Lambda=integrate(integrand, lower = 0, upper = t)$value
      Lambda+log(1-U[i])
    }
    r1 <- suppressWarnings(try(uniroot(f,  lower = 0, upper = 4), silent=TRUE))
    if (class(r1) == "try-error"){    
      pre_time[i]=4
    }
    else pre_time[i]=uniroot(f,  lower = 0, upper = 4)$root
  }
  
  pre_censoring=runif(N,0,3)
  pre_censoring=pre_censoring*(pre_censoring<3)+3*(pre_censoring>=3)
  tcens=(pre_censoring<pre_time) # censoring indicator
  delta=1-tcens
  time=pre_time*(delta==1)+pre_censoring*(delta==0)
  
  mean(delta)  # 0.496
  
  delta = delta[order(time)]
  z = z[order(time),]
  time = time[order(time)]
  time2=time[delta==1]
  
  df = cbind(time, delta, z)
  
  write.table(df,paste('../data/cox_time_varying_model/bdd_cox_N', N, '_',loop,'.csv', sep=''),row.names=FALSE,col.names=FALSE, sep=",") 
  
  setTxtProgressBar(pb, loop / nloop)
}
