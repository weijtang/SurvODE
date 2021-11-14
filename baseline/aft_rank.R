# this is fitting rank-based estimators for the semi-parametric aft model through aftgee package.
args = commandArgs(trailingOnly=TRUE) # sample size

library(R.matlab)
library(survival)
library(aftgee)

N=as.integer(args[1])
setting = 3

est_beta_aft = NULL
est_se_aft = NULL
semi_runtime = NULL

pb = txtProgressBar(style = 3)
for (i in 1:1000){
  print(i)
  #### load data
  survdata = readMat(paste('../data/survode/simudata_N', N, '_seed',i,'_setting', setting,'.mat', sep=''))
  x = survdata$x
  time = survdata$time
  delta = survdata$delta
  # cox model
  df = data.frame(x, time, delta)
  
  ############ fit semi-parametric aft ########################
  old_time=proc.time()
  rk.srrISMB <- aftsrr(Surv(time, delta) ~ ., data = data.frame(df), se = "ISMB") # rank based
  old_time=proc.time()- old_time
  semi_runtime=c(semi_runtime, old_time[3]) # record runtime
  
  # estimated beta
  beta_i = - rk.srrISMB$beta
  se_i = sqrt(diag(rk.srrISMB$covmat$ISMB))
  est_beta_aft = rbind(est_beta_aft, beta_i)
  est_se_aft = rbind(est_se_aft, se_i)
  
  setTxtProgressBar(pb, i / 1000)
}
save(semi_runtime, est_beta_aft, est_se_aft, file=paste('../res/baseline/bdd_aft_aftgee_N',N,'_setting', setting,'.RData', sep=''))

############### summary #################
setting = 3
Ns = c(1, 2, 4, 8)*1000
res_all = NULL
for (N in Ns){
  load(paste('../res/baseline/bdd_aft_aftgee_N',N,'_setting', setting,'.RData', sep=''))
  res_n= NULL
  temp = est_beta_aft
  res_n = colMeans(temp) - c(1,1,1)
  res_n = rbind(res_n,apply(temp, 2, sd))
  all_se = est_se_aft
  res_n = rbind(res_n,colMeans(all_se))
  low_ci = temp - 1.96*all_se
  up_ci = temp + 1.96*all_se
  true_beta = matrix(rep(c(1,1,1), 1000), ncol=3, byrow=TRUE)
  res_n = rbind(res_n,colMeans((low_ci < true_beta)*(up_ci > true_beta)))
  print(N)
  print(res_n)
  
  if (N==1000){
    runtime_base = mean(semi_runtime)
  }
  semi_runtime = semi_runtime / runtime_base
  temp = c(N, mean(semi_runtime),sd(semi_runtime))
  res_all = rbind(res_all, temp)
}
print(res_all)
write.table(res_all,paste('../res/baseline/bdd_aftgee_runtime.csv', sep=''),row.names=FALSE,col.names=FALSE, sep=",") 

