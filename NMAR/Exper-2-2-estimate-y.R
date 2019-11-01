##### Case d -compare all estimators -
rm(list=ls())
n <- 100000
ite <- 200

y_list <- seq(1,1,length.out = ite)

x_cate <- 2
for(iteite in 1:ite){
  print(iteite)
  x <- rmultinom(n,1,seq(1.0/x_cate,1.0/x_cate,length.out = x_cate))
  x1 <- rbinom(n,1,0.5)
  for(i in 1:n){
    x1[i] <- x[,i] %*% seq(0,x_cate-1.0,length.out = x_cate)
  }
  x2 <- runif(n,-1,1)
  y <- seq(0,0,length.out = n)
  delta <- seq(0,0,length.out = n)
  marge_pipi<- seq(0,0,length.out = n)
  mu <- -1.0-0.4*x1+0.5*x2*x2
  sigma_ <- 1.0
  alpha <- -0.3
  beta <-  -0.3
  beta2 <- -0.2
  gamma <- 0.5
  true_gamma <- gamma
  
  
  ### Marge PS
  phi <- sigma_*sigma_
  aaa <- (0.5*(gamma*phi+mu)*(gamma*phi+mu)-0.5*mu*mu)/phi
  for(i in 1:n){
    marge_pipi[i] <- 1.0/(1.0+exp(alpha+beta*x1[i]+beta2*x[i]*x[i]+aaa[i]))
    delta[i] <- rbinom(1,1,marge_pipi[i])
  }
  #### PS model
  sigmoid <- function(xxx){
    return(1.0/(1.0+exp(alpha+beta*x[i]+beta2*x[i]*x[i]+gamma*y[i])))
  }
  
  ### 
  for(i in 1:n){
    if(delta[i]==1){
      y[i] <- mu[i]+ rnorm(1,0,sigma_)
    }
    else{
      y[i] <- mu[i]+phi*gamma+rnorm(1,0,sigma_)
    }
  }
  
  likelihood_outcome <- function(gamma.){
    aaa <- 1.0/sqrt(gamma.[1]*gamma.[1])*exp(-0.5*(y[delta==1]-(gamma.[2]+gamma.[3]*x1[delta==1]+gamma.[4]*x2[delta==1]*x2[delta==1]))**2/(gamma.[1]*gamma.[1]))
    return(-sum(log(aaa)))
  }
  boku2 <- optim(c(1.0,-1.7,0.4,0.5),likelihood_outcome)
  mu <- (boku2$par[2]+boku2$par[3]*x1+boku2$par[4]*x2*x2)
  sigma_ <- boku2$par[1]
  
  y_list[iteite] <- mean(y)
}
print(mean(y_list))