###rm(list=ls())
library(gdata)
library(np)
### Mechanism 1
ite <- 100
y_list <- seq(1,1,length.out = ite) 
y_list2 <- seq(1,1,length.out = ite) 

for(iteite in 1:ite){
  print("ite")
  print(iteite)
  aaa <- read.xls("data1.xls",header=TRUE)
  aaa1 <- na.omit(aaa)
  nn <- 2506
  true_mean <- mean(aaa1$average.monthly.income.this.year..Korean.Won.10.000.)
  y <- aaa1$average.monthly.income.this.year..Korean.Won.10.000.
  y <- y/100.0
  x2 <- seq(1,1,length.out = nn)
  for(i in 1:nn){
    if(aaa1$age[i] < 35){
      x2[i] <- 0
    }
    else if (aaa1$age[i] < 51){
      x2[i] <- 1
    }
    else{
      x2[i] <- 2
    }
  }
  x1 <- aaa1$aveaage.monthly.income.in.the.previous.year......Korean.Won.10.000./100
  x3 <- aaa1$gender
  x4 <- aaa1$level.of.education
  
  #### PS model
  sigmoid <- function(xxx,yyy){
    return(1.0/(1.0+exp(-1.2-0.5*xxx-0.00*xxx*xxx+0.6*yyy)))
  }
  
  ### 
  delta <- seq(0,0,length.out = nn)
  for(i in 1:nn){
    delta[i] <- rbinom(1,1,sigmoid(x1[i],y[i]))
  }
  print(sum(delta==1)/2506)
  
  #### True_mean
  true_mean <- mean(y)
  ###y[delta==0] <- 0.0
  fit<- lm(y[delta==1]~x1[delta==1]*x2[delta==1]+x1[delta==1]*x3[delta==1]+x1[delta==1]*x4[delta==1]+x2[delta==1]*x3[delta==1]+x2[delta==1]*x4[delta==1]+x3[delta==1]*x4[delta==1], data=aaa1)
  AIC(fit)
  # ### Three interaction
  # fit<- lm(y[delta==1]~x1[delta==1]*x2[delta==1]+x1[delta==1]*x3[delta==1]+x1[delta==1]*x4[delta==1]+x2[delta==1]*x3[delta==1]+x2[delta==1]*x4[delta==1], data=aaa1)
  # AIC(fit)
  # fit<- lm(y[delta==1]~x1[delta==1]*x2[delta==1]+x1[delta==1]*x3[delta==1]+x1[delta==1]*x4[delta==1]+x2[delta==1]*x3[delta==1]+x3[delta==1]*x4[delta==1], data=aaa1)
  # AIC(fit)
  # fit<- lm(y[delta==1]~x1[delta==1]*x2[delta==1]+x1[delta==1]*x3[delta==1]+x1[delta==1]*x4[delta==1]+x2[delta==1]*x4[delta==1]+x3[delta==1]*x4[delta==1], data=aaa1)
  # AIC(fit)
  # fit<- lm(y[delta==1]~x1[delta==1]*x2[delta==1]+x1[delta==1]*x3[delta==1]+x2[delta==1]*x3[delta==1]+x2[delta==1]*x4[delta==1]+x3[delta==1]*x4[delta==1], data=aaa1)
  # AIC(fit)
  # fit<- lm(y[delta==1]~x1[delta==1]*x2[delta==1]+x1[delta==1]*x4[delta==1]+x2[delta==1]*x3[delta==1]+x2[delta==1]*x4[delta==1]+x3[delta==1]*x4[delta==1], data=aaa1)
  # AIC(fit)
  # fit<- lm(y[delta==1]~x1[delta==1]*x3[delta==1]+x1[delta==1]*x4[delta==1]+x2[delta==1]*x3[delta==1]+x2[delta==1]*x4[delta==1]+x3[delta==1]*x4[delta==1], data=aaa1)
  # #### No interaction
  # AIC(fit)
  # fit<- lm(y[delta==1]~x1[delta==1]+x2[delta==1]+x3[delta==1]+x4[delta==1], data=aaa1)
  
  ### Missing rate
  
  tttt <- fit$coefficients
  mu <- tttt[1]+tttt[2]*x1+tttt[3]*x2+tttt[4]*x3+tttt[5]*x4+tttt[6]*x1*x2+tttt[7]*x1*x3+tttt[8]*x1*x4+tttt[9]*x2*x3+tttt[10]*x2*x4+tttt[11]*x3*x4
  sigma_ <- summary(fit)$sigma
  
  gamma <- 0.6
  ccc10 <- npksum(txdat = x1, tydat = delta*exp(gamma*y),ckerorder=2,bws=0.25)
  ccc11 <-  npksum(txdat = x1, tydat = 1-delta,ckerorder=2,bws=0.25)
  pseudo_pi <- as.vector(ccc10$ksum/(ccc10$ksum+ccc11$ksum*exp(gamma*y)))
  plot(x1, 1.3+0.2*x1+0.3*x1**0.5,ylim = c(-5,10))
  par(new=TRUE)
  plot(x1,log(as.vector(ccc10$ksum/ccc11$ksum)),col="green",ylim = c(-5,10))
  
  likelihood <- function(gamma.,y.=y){
    pseudo_pi <- 1.0/(1.0+exp(gamma.[1]+gamma.[2]*x1+gamma.[3]*(x1*x1)+gamma.[4]*y.))
    aaa <- delta*log(pseudo_pi)+(1-delta)*log(1-pseudo_pi)
    return(-sum(aaa))
  }
  boku2 <- optim(c(-1.3,-0.2,0.3,0.3),likelihood)
  #####print(boku2$par)
  gamma_estimate <- boku2$par[4]
  print(gamma_estimate)
  # #### Kott
  # 
  # likelihood2 <- function(gamma.,y.=y){
  #   ###gamma.[1] <- alpha
  #   ####gamma.[2] <- beta
  #   pseudo_pi <- 1.0/(1.0+exp(gamma.[1]+gamma.[2]*x1+gamma.[3]*y.))
  #   ####pseudo_pi <- 1.0/(1.0+exp(-1.7-0.4*x1+gamma.*y.))
  #   aaa <- (delta/pseudo_pi-1.0)*x2
  #   aaa1 <- (delta/pseudo_pi-1.0)*x1
  #   aaa2 <- (delta/pseudo_pi-1.0)
  #   return(sum(aaa)*sum(aaa)+sum(aaa1)*sum(aaa1)+sum(aaa2)*sum(aaa2))
  # }
  # boku2 <- optim(c(-1.7,-0.4,0.5),likelihood2)
  # ####print(boku2$par)
  # gamma_estimate2 <- boku2$par[3]
  # 
  #### Optimal IV
  
  likelihood3 <- function(gamma.,y.=y){
    pseudo_pi <- 1.0/(1.0+exp(gamma.[1]+gamma.[2]*x1+gamma.[3]*y.))
    pseudo_pi2 <- 1.0/(1.0+exp(gamma.[1]+gamma.[2]*x1+gamma.[3]*(1.5*gamma.[3]*sigma_*sigma_+mu)))
    aaa <- (delta/pseudo_pi-1.0)*pseudo_pi2*(gamma.[3]*sigma_*sigma_+mu)
    aaa1 <- (delta/pseudo_pi-1.0)*pseudo_pi2
    aaa2 <- (delta/pseudo_pi-1.0)*pseudo_pi2*(x1)
    return(sum(aaa)*sum(aaa)+sum(aaa1)*sum(aaa1)+sum(aaa2)*sum(aaa2))
  }
  boku3 <- optim(c(-2.0,-0.3,0.3),likelihood3)
  gamma_estimate3 <- boku3$par[3]
  ######print(gamma_estimate3)
  
  #### Proposed
  
  likelihood4 <- function(gamma.){
    ccc1 <- npksum(txdat = x1, tydat = delta*exp(gamma.*y),ckerorder=2,bws=0.25)
    ccc2 <-  npksum(txdat = x1, tydat = 1-delta,ckerorder=2,bws=0.25)
    pseudo_pi <- as.vector(ccc1$ksum/(ccc1$ksum+ccc2$ksum*exp(gamma.*y)))
    pseudo_pi2 <- as.vector(ccc1$ksum/(ccc1$ksum+ccc2$ksum*exp(gamma.*(1.5*gamma.*sigma_*sigma_+mu))))
    aaa <- (delta/pseudo_pi-1.0)*pseudo_pi2*(gamma.*sigma_*sigma_+mu)
    return(sum(aaa)*sum(aaa))
  }
  
  boku4 <- optimise(likelihood4,c(-2.0,2.0))
  gamma_estimate4 <- boku4$minimum
  print("estimate gamma")
  print(boku4$minimum)
  
  #### IPW 
  
  ccc1 <- npksum(txdat = x1, tydat = delta*exp(boku4$minimum*y),bws = 0.3)
  ccc2 <- npksum(txdat = x1, tydat = 1-delta,bws = 0.3)
  g_list2 <- as.vector(ccc1$ksum/ccc2$ksum)
  pseudo_pi <- as.vector(g_list2/(g_list2+exp(boku4$minimum*y)))
  y_estimate4_ipw<- mean(na.omit(y/pseudo_pi))
  print(y_estimate4_ipw)
  
  ### DB
  
  ccc1 <- npksum(txdat = x1, tydat = delta*exp(boku4$minimum*y),bws = 0.3)
  ccc2 <- npksum(txdat = x1, tydat = 1-delta,bws = 0.3)
  g_list2 <- as.vector(ccc1$ksum/ccc2$ksum)
  pseudo_pi <- as.vector(g_list2/(g_list2+exp(boku4$minimum*y)))
  y_estimate4_db <- mean(na.omit((y-mu)*delta/pseudo_pi+mu))
  print(y_estimate4_db)
  
  ### outcome
  
  
  ## Wang&Shao
  likelihood5 <- function(gamma.){
    ccc1 <- npksum(txdat = x1, tydat = delta*exp(gamma.*y),ckerorder=2,bws=0.25)
    ccc2 <-  npksum(txdat = x1, tydat = 1-delta,ckerorder=2,bws=0.25)
    pseudo_pi <- as.vector(ccc1$ksum/(ccc1$ksum+ccc2$ksum*exp(gamma.*y)))
    aaa <- (delta/pseudo_pi-1.0)*x2
    return(sum(aaa)*sum(aaa))
  }
  boku5 <- optimise(likelihood5,c(-2.0,2.0))
  gamma_estimate5 <- boku5$minimum
  print("estimate gamma")
  print(boku5$minimum)
  
  
  #### IPW 
  ccc1 <- npksum(txdat = x1, tydat = delta*exp(boku5$minimum*y),bws = 0.3)
  ccc2 <- npksum(txdat = x1, tydat = 1-delta,bws = 0.3)
  g_list2 <- as.vector(ccc1$ksum/ccc2$ksum)
  pseudo_pi <- as.vector(g_list2/(g_list2+exp(boku5$minimum*y)))
  y_estimate5_ipw<- mean(na.omit(y/pseudo_pi))
  print(y_estimate5_ipw)
  
  ### DB
  
  ccc1 <- npksum(txdat = x1, tydat = delta*exp(boku5$minimum*y),bws = 0.3)
  ccc2 <- npksum(txdat = x1, tydat = 1-delta,bws = 0.3)
  g_list2 <- as.vector(ccc1$ksum/ccc2$ksum)
  pseudo_pi <- as.vector(g_list2/(g_list2+exp(boku5$minimum*y)))
  y_estimate5_db <- mean(na.omit((y-mu)*delta/pseudo_pi+mu))
  print(y_estimate5_db)
  
  #### Variance estimator
  ccc1 <- npksum(txdat = x1, tydat = delta*exp(gamma*y),bws = 0.3)
  ccc2 <- npksum(txdat = x1, tydat = 1-delta,bws = 0.3)
  g_list2 <- as.vector(ccc1$ksum/ccc2$ksum)
  pseudo_pi <- as.vector(g_list2/(g_list2+exp(gamma*y)))
  
  fit2 <- lm(y[delta==1]~x1[delta==1], data=aaa1)
  tttt2 <- fit2$coefficients
  mu2 <- tttt2[1]+tttt2[2]*x1
  A2 <- (1-pseudo_pi)*(y-mu2)*(mu-mu2)
  A2 <- A2[A2<5]
  B2 <- (1-pseudo_pi)/pseudo_pi*(mu2-mu)*(mu2-mu)
  B2 <- B2[B2<5]
  C2 <- (1-pseudo_pi)*(y-mu)*(y-mu)
  C2 <- C2[C2<5]
  A2_mean <- mean(na.omit(A2))
  B2_mean <- mean(na.omit(B2))
  C2_mean <- mean(na.omit(C2))
  
  print(((var(na.omit((y-mu)*delta/pseudo_pi+mu))+C2_mean*B2_mean/A2_mean)/nn)**0.5)
  y_list[iteite] <- y_estimate4_db
  y_list2[iteite] <- y_estimate5_db
}
print(sqrt(var(y_list)))
print(sqrt(var(y_list2)))

