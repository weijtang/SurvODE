# this is fitting mple for the semi-parametric cox model.
args = commandArgs(trailingOnly=TRUE)

library(R.matlab)
library(survival)

N=as.integer(args[1])
a = seq(0, 1.5, length.out=101)

NR.time=NULL
est_beta_cox = NULL
est_se_cox = NULL
bh = NULL

setting = 1

for (i in 1:1000){
  print(i)
  survdata = readMat(paste('../data/survode/simudata_N', N, '_seed',i,'_setting', setting,'.mat', sep=''))
  x = survdata$x
  time = survdata$time
  delta = survdata$delta
  # cox model
  df = data.frame(x, time, delta)
  old_time3=proc.time()
  cox_fit = coxph(Surv(time, delta) ~ ., data = df)
  old_time3=proc.time()- old_time3
  NR.time=c(NR.time, old_time3[3])
  
  beta_i = cox_fit$coefficients
  var_i = cox_fit$var
  se_i = sqrt(diag(var_i))
  est_beta_cox = rbind(est_beta_cox, beta_i)
  est_se_cox = rbind(est_se_cox, se_i)
  
  temp = basehaz(cox_fit,centered=FALSE)
  px = temp$time
  idx = findInterval(a, px)
  idx[which(idx==0)] = 1
  length(idx)
  est_bh = temp$hazard
  bh = rbind(bh, t(est_bh[idx]))
  
}

save(NR.time, est_beta_cox, est_se_cox, bh,file=paste('../res/baseline/bdd_cox_mple_N',N,'_setting', setting,'.RData', sep=''))

############### summary #################
setting = 1
Ns = c(1, 2, 4, 8)*1000
for (N in Ns){
  load(paste('../res/baseline/bdd_cox_mple_N',N,'_setting', setting,'.RData', sep=''))
  res_n= NULL
  temp = est_beta_cox
  res_n = colMeans(temp) - c(1,1,1)
  res_n = rbind(res_n,apply(temp, 2, sd))
  all_se = est_se_cox
  res_n = rbind(res_n,colMeans(all_se))
  low_ci = temp - 1.96*all_se
  up_ci = temp + 1.96*all_se
  true_beta = matrix(rep(c(1,1,1), 1000), ncol=3, byrow=TRUE)
  res_n = rbind(res_n,colMeans((low_ci < true_beta)*(up_ci > true_beta)))
  print(N)
  print(res_n)
}
