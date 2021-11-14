library(readr)
N = 1000
NR.time = read_csv(paste('../res/baseline/bdd_cox_time_varying_runtime_N',N,'.csv', sep=''), col_names = FALSE) 
all_alpha = read_csv(paste('../res/baseline/bdd_cox_time_varying_est_tv_N',N,'.csv', sep=''), col_names = FALSE)
all_beta = read_csv(paste('../res/baseline/bdd_cox_time_varying_est_tinv_N',N,'.csv', sep=''), col_names = FALSE) 
all_se = read_csv(paste('../res/baseline/bdd_cox_time_varying_se_tinv_N',N,'.csv', sep=''), col_names = FALSE) 
NR.time = as.matrix(NR.time)
all_alpha = as.matrix(all_alpha)
all_beta = as.matrix(all_beta)
all_se = as.matrix(all_se)
# print computation time
mean(NR.time)
sd(NR.time)


res = NULL
# print estimated beta
round(colMeans(all_beta) - c(1,-1,-1,1), 3)
round(apply(all_beta, 2, sd), 3)
true_beta = matrix(rep(c(1,-1,-1,1), 1000), ncol=4, byrow=TRUE)

# print ese and cp
round(colMeans(all_se),3)
low_ci = all_beta - 1.96*all_se
up_ci = all_beta + 1.96*all_se
round(colMeans((low_ci < true_beta)*(up_ci > true_beta)), 3)

# print IMSE of time varying effects
px = seq(0, 2.0, length.out=100)
true_alpha = sin(3*pi*px/4)
imse = NULL
for (i in 1:nrow(all_alpha)){
  imse[i] = sum((all_alpha[i,] - true_alpha)^2) * 2/99
}
mean(imse)
sd(imse)