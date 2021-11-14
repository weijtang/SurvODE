# this is fitting the cox model through flexsurv package.
args = commandArgs(trailingOnly=TRUE) # sample size

library(R.matlab)
library(survival)
library(splines)
library(flexsurv)

N=as.integer(args[1])
a = seq(0, 1.5, length.out=101)

NR.time=NULL
all_cumhaz = NULL
all_haz=NULL
all_beta=NULL
all_se= NULL

setting = 1
pb = txtProgressBar(style = 3)
for (i in 1:1000){
  print(i)
  survdata = readMat(paste('../data/survode/simudata_N', N, '_seed',i,'_setting', setting,'.mat', sep=''))
  x = survdata$x
  time = survdata$time
  delta = survdata$delta
  # cox model
  df = data.frame(x, time, delta)
  num_k = ceiling(length(unique(time))^(1/5)) - 1
  old_time3=proc.time()
  cox_fit = flexsurvspline(Surv(time, delta) ~ X1+X2+X3, data = df, k=num_k, scale="hazard", timescale="log")
  old_time3=proc.time()- old_time3
  NR.time=c(NR.time, old_time3[3])
  
  # estimated cumulative hazard
  px = seq(0.01, 1.5, length.out=100)
  est_cumhaz = predict(cox_fit,newdata=data.frame(X1=0, X2=0, X3=0), type="cumhaz", times=px)$.pred[[1]]
  all_cumhaz = rbind(all_cumhaz, t(est_cumhaz[,2]))
  est_haz = predict(cox_fit,newdata=data.frame(X1=0, X2=0, X3=0), type="hazard", times=px)$.pred[[1]]
  all_haz = rbind(all_haz, t(est_haz[,2]))
  
  # estimated beta
  all_beta = rbind(all_beta, cox_fit$coefficients[(num_k+3):length(cox_fit$coefficients)])
  all_se = rbind(all_se, sqrt(diag(cox_fit$cov)[(num_k+3):length(cox_fit$coefficients)]))
  setTxtProgressBar(pb, i / 1000)
}
save(NR.time, all_beta, all_se, all_haz, all_cumhaz, file=paste('../res/baseline/bdd_cox_flexsurv_N',N,'_setting', setting,'.RData', sep=''))

############### summary #################
setting = 1
Ns = c(1, 2, 4, 8)*1000
res_all = NULL
for (N in Ns){
  load(paste('../res/baseline/bdd_cox_flexsurv_N',N,'_setting', setting,'.RData', sep=''))
  temp = all_beta
  res_n = colMeans(temp) - c(1,1,1)
  res_n = rbind(res_n,apply(temp, 2, sd))
  res_n = rbind(res_n,colMeans(all_se))
  low_ci = temp - 1.96*all_se
  up_ci = temp + 1.96*all_se
  true_beta = matrix(rep(c(1,1,1), 1000), ncol=3, byrow=TRUE)
  res_n = rbind(res_n,colMeans((low_ci < true_beta)*(up_ci > true_beta)))
  print(N)
  print(res_n)
  
  # print IMSE of cumulative hazard
  px = seq(0.01, 1.5, length.out=100)
  true_cumhaz = px^4/4
  imse_cumhaz = NULL
  for (i in 1:nrow(all_cumhaz)){
    imse_cumhaz[i] = sum((all_cumhaz[i,] - true_cumhaz)^2) * 1.5/99
  }
  # print IMSE of hazard
  px = seq(0.01, 1.5, length.out=100)
  true_haz = px^3
  imse_haz = NULL
  for (i in 1:nrow(all_haz)){
    imse_haz[i] = sum((all_haz[i,] - true_haz)^2) * 1.5/99
  }
  imse_haz = imse_haz 
  imse_cumhaz = imse_cumhaz 
  if (N==1000){
    runtime_base = NR.time
  }
  NR.time = NR.time / runtime_base
  temp = c(N, mean(imse_haz), sd(imse_haz),
           mean(imse_cumhaz), sd(imse_cumhaz),
           mean(NR.time),sd(NR.time))
  res_all = rbind(res_all, temp)
}
print(res_all)
write.table(res_all,paste('../res/baseline/bdd_cox_flexsurv_all.csv', sep=''),row.names=FALSE,col.names=FALSE, sep=",") 
