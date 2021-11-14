args = commandArgs(trailingOnly=TRUE) # sample size

library(survival)
library(mvtnorm)
library(splines2)
p=5
N=as.integer(args[1])
tol=10^(-6)

df=3

NR.time=NULL
all_alpha=NULL
all_beta=NULL
all_se= NULL



loop=0
pb = txtProgressBar(style = 3)
nloop=1000
while (loop<nloop) {
  loop=loop+1
  set.seed(loop)
  ######### load data #############
  df = read.table(paste('../data/cox_time_varying_model/bdd_cox_N', N, '_',loop,'.csv', sep=''),row.names=NULL,sep=",")
  time = df$V1
  delta = df$V2
  z = as.matrix(df[,3:(2+p)])
  time2=time[delta==1]
  
  #initial values
  vfit <- survival::coxph(Surv(time,delta) ~z)
  constant_beta<-vfit$coef
  
  knot = ceiling(length(unique(time))^(1/5)) + 3
  ZZ=z[,2]
  
  old_time3=proc.time()
  output= coxph(Surv(time, delta)~ tt(ZZ)+z[,-2], 
                tt=function(ZZ,time,...) {mtrx = ZZ
                t=time
                knots=quantile(time2,prob=seq(1:(knot-4))/(knot-3),type=2)
                X=bSpline(t,knots=knots,intercept=TRUE,degree=3L,Boundary.knots=c(0, max(time))) 
                ZZ * X
                })
  old_time3=proc.time()- old_time3
  print("done")
  
  NR.time=c(NR.time, old_time3[3])
  
  # estimated time varying effects
  px = seq(0, 2.0, length.out=100)
  knots = quantile(time2,prob=seq(1:(knot-4))/(knot-3),type=2)
  est_alpha = bSpline(px,knots=knots,intercept=TRUE,degree=3L,Boundary.knots=c(0, max(time)))%*%(as.matrix(output$coefficients[1:(length(knots)+4)]))
  all_alpha = rbind(all_alpha, t(est_alpha))
  
  # estimated beta
  all_beta = rbind(all_beta, output$coefficients[(length(knots)+5):length(output$coefficients)])
  all_se = rbind(all_se, sqrt(diag(output$var)[(length(knots)+5):length(output$coefficients)]))
  setTxtProgressBar(pb, loop / nloop)
}

write.table(NR.time,paste('../res/baseline/bdd_cox_time_varying_runtime_N',N,'.csv', sep=''),row.names=FALSE,col.names=FALSE, sep=",") 
write.table(all_alpha,paste('../res/baseline/bdd_cox_time_varying_est_tv_N',N,'.csv', sep=''),row.names=FALSE,col.names=FALSE, sep=",") 
write.table(all_beta,paste('../res/baseline/bdd_cox_time_varying_est_tinv_N',N,'.csv', sep=''),row.names=FALSE,col.names=FALSE, sep=",") 
write.table(all_se,paste('../res/baseline/bdd_cox_time_varying_se_tinv_N',N,'.csv', sep=''),row.names=FALSE,col.names=FALSE, sep=",") 

# print computation time
mean(NR.time)
sd(NR.time)

# print estimated beta
colMeans(all_beta) - c(1, -1, -1, 1)
apply(all_beta, 2, sd)
colMeans(all_se)

# print IMSE of time varying effects
px = seq(0, 2.0, length.out=100)
true_alpha = sin(3*pi*px/4)
imse = NULL
for (i in 1:nrow(all_alpha)){
  imse[i] = mean((all_alpha[i,] - true_alpha)^2) * 3
}
mean(imse)
sd(imse)

pdf(file=paste('../res/baseline/bdd_time_varying_effects_N',N, '.pdf', sep=''))
plot(px,colMeans(all_alpha), col="red", ylim=c(-2,2), xlab="Time", ylab=expression(beta[1](t)), main="", lty=1, cex=1, type='l', cex.lab=1, cex.axis=1,lwd=2)
lines(px, true_alpha, lty=2, col="blue")

alpha_low=rep(0, ncol(all_alpha))
alpha_up=rep(0, ncol(all_alpha))         
for (i in 1:ncol(all_alpha))
{
  alpha_low[i]=quantile(all_alpha[,i], probs = 0.025)
  alpha_up[i]=quantile(all_alpha[,i], probs = 0.975)
}
lines(px, alpha_low, lty=5, col="orange", lwd=2)
lines(px, alpha_up, lty=5, col="orange", lwd=2)
legend ("bottomleft", legend=c('Estimate'
                               , 
                               'True value'
                               ,
                               '95% Percentile'), yjust = 1,text.width = strwidth("1,000,000,000"),lty=c(1,2, 5), pt.cex = 1, cex=0.9,col=c("black", "blue", "orange"),lwd=2)
dev.off()


